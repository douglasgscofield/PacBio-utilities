#!/usr/bin/env perl

use strict;
use warnings;
# base modules
use POSIX;
use Getopt::Long;
use List::Util qw/sum/;
# BioPerl modules
use Bio::Seq;
use Bio::SeqIO;

my $o_mappingquality = 0;
my $o_mincov = 0;
my $o_minqual = 0;
my $o_minqualcov = 0;
my $o_variantsonly = 0;
my $o_printbasesquals = 0;
my $o_noposcheck = 0;
my $o_help;
sub print_usage_and_exit($$);
sub merge_pileup_columns($);
sub indel_targets();
sub indel_apply();

################
# main loop

my %command_values = ( "indel-targets" => "indel-targets",
                       "indel-apply" => "indel-apply");
my $command = shift @ARGV;

if (not defined $command or not defined $command_values{$command}) {
    print STDERR "$0 indel-targets ...\n";
    print STDERR "$0 indel-apply ...\n";
    exit 0;
} elsif ($command eq "indel-targets") {
    indel_targets();
} elsif ($command eq "indel-apply") {
    indel_apply();
}

################

sub print_usage_and_exit($$) {
    my $usage = shift;
    my $msg = shift;
    print "$msg\n" if $msg;
    print $usage;
    exit 0;
}

################

my $do_merge = 0;
sub merge_pileup_columns($) {
    # each BAM has 4 columns: coverage, bases, quality, mapping quality
    my $l = shift;
    my ($cvg, $base, $qual, $mapqual) = ( 0, "", "", "" );
    my ($n_cols, $n_bams) = ( scalar(@$l), 0 );
    my $i = 3; # start at 4th column (coverage)
    while ($i < $n_cols) {
        if ($l->[$i] > 0) {  # coverage > 0, 4 columns
            $cvg     += $l->[$i];
            $base    .= $l->[$i + 1];
            $qual    .= $l->[$i + 2];
            $mapqual .= $l->[$i + 3];
            $i += 4;
        } else {  # no coverage here, 3 columns
            $i += 3;
        }
        ++$n_bams;
    }
    if (not $do_merge) {
        print STDERR "Merging columns from $n_bams BAM files ...\n";
        $do_merge = 1;
    }
    return ($l->[0], $l->[1], $l->[2], $cvg, $base, $qual, $mapqual);
}

################

sub indel_targets() {
    my $o_fasta;
    my $o_stdin;
    my $o_indelfrac = 0.80;
    my $o_indelmincov = 5;
    my $o_includebadindels = 0;
    my $o_maxindelsize = 1000;
    my $o_basequaloffset = 33; # by default, assume Phred+33 quality scores
    my $o_help = 0;
    my $o_samtools;
    my $o_allreads;
    my $min_mapqual = 1;
    my $min_basequal = 13;
    my $samtools_filterflags = "--ff ".(256 + 1024);  # not primary, is duplicate
    my $o_allbases;
    my $samtools = "samtools";


    my $usage = "
pacbio-util.pl indel-targets -f FASTA FILE1.bam [ FILE2.bam ... ]

    -f FILE, --fasta FILE    PacBio assembly in Fasta format

    FILE1.bam [ FILE2.bam ]  BAM files of reads mapped to

    -                        Read pileup from stdin rather than running samtools

    --indel-frac FLOAT       do not report indels at a position if the fraction of 
                             reads containing them is below FLOAT [default $o_indelfrac]
    --indel-min-cov INT      Minimum read coverage to consider an indel [default $o_indelmincov]
    --include-bad-indels     Include indels in target set that do not pass our criteria,
                             these lines are marked with 'bad' in the sixth column
    --max-indel-size INT     Maximum indel size to accept, beyond this the position is
                             considered in error [default $o_maxindelsize]
    --base-qual-offset INT   Offset of base quality ASCII value from 0 quality [default $o_basequaloffset]
    --samtools FILE          Location of samtools executable [default $samtools]
    --all-reads              Include all reads.  Default is to exclude multiply-mapped
                             reads (mapQ 0), alignments that are not the primary (SAM
                             flag 256), and duplicate alignments (flag 1024).
    --all-bases              Include all bases.  Default is to exclude bases below
                             quality 13.

    --help, -?               help message

";

    print_usage_and_exit($usage, "") if not @ARGV;
    GetOptions(""                   => \$o_stdin,
               "fasta=s"            => \$o_fasta,
               "indel-frac=f"       => \$o_indelfrac,
               "indel-min-cov=i"    => \$o_indelmincov,
               "include-bad-indels" => \$o_includebadindels,
               "max-indel-size=i"   => \$o_maxindelsize,
               "base-qual-offset=i" => \$o_basequaloffset,
               "samtools=s"         => \$o_samtools,
               "all-reads"          => \$o_allreads,
               "all-bases"          => \$o_allbases,
               "help|?"             => \$o_help
    ) or print_usage_and_exit($usage, "unknown option");
    print_usage_and_exit($usage, "") if $o_help;
    my @BAM = grep { -f or die "cannot find BAM: $_" } @ARGV;
    my $n_BAM = scalar(@BAM);

    # find samtools
    if ($o_samtools) {
        die "cannot execute '$o_samtools'" if not -x $o_samtools;
        $samtools = $o_samtools;
    } else {
        my $t = qx/$samtools 2>&1/;
        die "cannot execute '$samtools'" if not defined $t;
    }
    # set up samtools argument list
    my $samtools_args = "-s -B -d 10000 -L 10000";
    $samtools_args .= " $samtools_filterflags" if not $o_allreads;
    $min_mapqual = 0 if $o_allreads;
    $samtools_args .= " -q $min_mapqual";
    $min_basequal = 0 if $o_allbases;
    $samtools_args .= " -Q $min_basequal";
    # open samtools on a pipe
    my $samtools_pipe = "$samtools mpileup $samtools_args -f $o_fasta ".join(" ", @BAM)." |";
    print STDERR "pipe: $samtools_pipe\n";
    open(PILEUP, $samtools_pipe) or die "Could not initiate samtools pipe: $!";

    print STDOUT "#assembly:$o_fasta\n";

    while (<PILEUP>) {
        next if m/^#/;
        chomp;
        my @l = split /\t/;
        my ($ref, $pos, $refcall, $cvg, $base, $qual, $mapqual ) = merge_pileup_columns(\@l);
        die "Coverage contains non-numeric values: $cvg" if not isdigit($cvg);
        if ($refcall eq "*") {
            print STDERR "skipping reference deletion '*', should this happen??\n";
            next; # skip indel lines
        }
        my %indels;
        my %indel_ops;
        my %indel_lengths;
        my $in_reads = 0;
        my $del_reads = 0;
        my $mincvg_not_met = 1 if $cvg < $o_mincov;
        my $minqualcvg_not_met = 0;
        my $n_bases;
        if ($base =~ m/[\$\^\+-]/) {
            $base =~ s/\^.//g; #removing the start of the read segement mark
            $base =~ s/\$//g; #removing end of the read segment mark
            while ($base =~ m/([\+-]){1}(\d+)/g) {
                my $indel_type = $1;
                my $indel_len = $2;
                if ($indel_len > $o_maxindelsize) { # problem with pileup line?
                    print STDERR "line $. indel too big, $indel_len\n";
                    next;
                }
                $base =~ s/([\+-]{1}$indel_len.{$indel_len})//; # remove indel info from read base field
                my $indel_contents = $1;
                ++$indels{"indel, $indel_type, $indel_len bases, $indel_contents"};
                ++$indel_lengths{($indel_type eq "+" ? +$indel_len : -$indel_len)};
                ++$indel_ops{uc($indel_contents)};
                ++$in_reads if $indel_type eq "+";
                ++$del_reads if $indel_type eq "-";
            }
            $n_bases = length($base);
            die "n bases $n_bases != number of base qualities" if $n_bases != length($qual);
            # below will not be correct when we are adding a new read ^I, mapping qualities wil
            # contain the quality of the new read but the base and base qual fields will not
            #die "n bases $n_bases != number of map qualities" if $n_bases != length($mapqual);
        }
        my $indel_reads = $in_reads + $del_reads;
        next if $indel_reads == 0;  # no indel
        my $indel_frac_pos = $indel_reads / $n_bases;
        my $indel_qual = ($indel_frac_pos >= $o_indelfrac and 
                          $indel_reads >= $o_indelmincov and
                          scalar(keys %indel_lengths) == 1 and
                          scalar(keys %indel_ops) == 1) ? "good" : "bad";

        next if not $o_includebadindels and $indel_qual eq "bad";
        my @out;
        push @out, ($ref, $pos, $refcall, $cvg, "indel-target");
        push @out, $indel_qual;
        push @out, (sprintf "%.4f", $indel_frac_pos);
        push @out, join(",", keys %indel_lengths);
        push @out, join(",", keys %indel_ops);
        push @out, $in_reads + $del_reads;
    	print STDOUT join("\t", @out), "\n";
    }
    close PILEUP;
}

################

sub find_ref_fasta_seq($);

my $ref_fasta_io;   # SeqIO object for reading assembly
my $ref_fasta_seq;  # current Seq object from assembly


sub indel_apply() {
    print STDERR "you've reached indel_apply\n";
    my $usage = "
pacbio-util.pl indel-apply -f FASTA -t TARGETS

Correct indels within FASTA by applying TARGETS.  The order of sequences in FASTA
and targets in TARGETS must be the same.

    -f | --fasta FILE        PacBio assembly in Fasta format, to be corrected.
                             Must be the same assembly used for
                             '$0 indel-targets -f FILE ...'
    -t | --targets TARGETS   List of indel targets to apply, output from
                             '$0 indel-targets -f FILE ...'

    --help, -?               help message

";

    # options for ...

    my $o_fasta;
    my $o_targets;
    my $o_stdin;

    print_usage_and_exit($usage, "") if not @ARGV;
    GetOptions(""          => \$o_stdin,
               "fasta=s"   => \$o_fasta,
               "targets=s" => \$o_targets
    ) or print_usage_and_exit($usage, "unknown option");

    open(TARGETS, "<$o_targets") or die "Could not open targets file '$o_targets': $!";
    my $target_assembly = <TARGETS>;  # first line is '#assembly:filename'
    die "First line of '$o_targets' does not contain assembly filename" if $target_assembly !~ m/^#assembly:/;
    chomp $target_assembly;
    $target_assembly = substr($target_assembly, 10);  # remove #assembly:
    if (not defined $o_fasta) {
        print STDERR "Opening target assembly '$target_assembly'\n" if not $o_fasta;
        $o_fasta = $target_assembly;
    } elsif ($o_fasta ne $target_assembly) {
        print STDERR "\n*** Warning: -f and target assembly do not match\n\n";
    }
    $ref_fasta_io = Bio::SeqIO->new(-file => "<$o_fasta", -format => "fasta") or 
                    die "Could not open fasta file '$o_fasta': $!";

    while (<TARGETS>) {
        next if m/^#/;
        chomp;
        # ?what is the target line format?
        # skip 'bad' indels
        # if the reference sequence is new, write out the previous one
        # get the indel to apply (generally a single-base insertion, yes?)
        # apply it to the reference sequence
    }
    # write out the reference sequence
}

sub find_ref_fasta_seq($) {
    # read fasta sequences until we load the requested one
    my $seq_name_to_find = shift;
    # if we return below, then we are already at the sequence of interest
    if (defined($ref_fasta_seq) and $ref_fasta_seq->id eq $seq_name_to_find) {
        #print STDERR "find_ref_fasta_seq: we are already at the reference sequence $seq_name_to_find\n";
        return;
    }
    while (my $s = $ref_fasta_io->next_seq()) {
        if ($s->id eq $seq_name_to_find) {
            $ref_fasta_seq = $s;
            #print STDERR "find_ref_fasta_seq: found new reference sequence $seq_name_to_find\n";
            return;
        }
    }
    die("couldn't find reference $seq_name_to_find");
}

