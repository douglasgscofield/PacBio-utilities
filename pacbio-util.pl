#!/usr/bin/env perl

use strict;
use warnings;
# base modules
use POSIX;
use File::Which;
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
my $o_basequaloffset = 33; # by default, assume Phred+33 quality scores
my $o_indels = 1;
my $o_indelmode = 1;
my $o_indelfrac = 0.80;
my $o_indeltoobig = 1000;

sub print_usage_and_exit($$);
sub merge_pileup_columns($);
sub indel_targets();
sub indel_apply();

################
# main loop

my $command = shift @ARGV;

if ($command eq "indel-targets") {
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

sub merge_pileup_columns($) {
    # each BAM has 4 columns: coverage, bases, quality, mapping quality
    my $l = shift;
    my ($cvg, $base, $qual, $mapqual) = ( 0, "", "", "" );
    my $n_cols = scalar(@$l);
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
    }
    return ($l->[0], $l->[1], $l->[2], $cvg, $base, $qual, $mapqual);
}

################

sub indel_targets() {
    my $usage = "
pacbio-util.pl indel-targets -f FASTA FILE1.bam [ FILE2.bam ... ]

    -f FILE, --fasta FILE    PacBio assembly in Fasta format
    FILE1.bam [ FILE2.bam ]  BAM files of reads mapped to

    --base-qual-offset INT   offset of base quality ASCII value from 0 quality [default $o_basequaloffset]

    --indel-frac FLOAT       do not report indels at a position if the fraction of 
                             reads containing them is below FLOAT [default $o_indelfrac]
    --indel-mode             track ONLY the presence of indels [default $o_indelmode]

    --help, -?               help message

";

    my $o_fasta;
    my $o_stdin;

    GetOptions(""          => \$o_stdin,
               "fasta=s"   => \$o_fasta
    ) or print_usage_and_exit($usage, "unknown option");
    my @BAM = grep { -f or die "cannot find BAM: $_" } @ARGV;
    my $n_BAM = scalar(@BAM);

    my $samtools = which('samtools') or die "$command: samtools not found";

    my $samtools_pipe = "$samtools mpileup -s -BQ0 -d 1000000 -L 1000000 -f $o_fasta ".join(" ", @BAM)." |";
    open(PILEUP, $samtools_pipe) or die "Could not initiate samtools pipe: $!";

    print STDOUT "#assembly:$o_fasta\n";

    while (<PILEUP>) {
        next if m/^#/;
        chomp;
        my @l = split /\t/;
        my ($ref, $pos, $refcall, $cvg, $base, $qual, $mapqual ) = merge_pileup_columns(\@l);
        die "Coverage contains non-numeric values: $cvg" if not isdigit($cvg);
        my $n_bases = length($base);
        die "number of bases bases does not match number of base qualities" if $n_bases != length($qual);
        die "number of bases bases does not match number of map qualities" if $n_bases != length($mapqual);
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
        my $note = "";
        if ($base =~ m/[\$\^\+-]/) {
            $base =~ s/\^.//g; #removing the start of the read segement mark
            $base =~ s/\$//g; #removing end of the read segment mark
            while ($base =~ m/([\+-]){1}(\d+)/g) {
                my $indel_type = $1;
                my $indel_len = $2;
                if ($indel_len > $o_indeltoobig) { # problem with pileup line?
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
        }
        next if $in_reads + $del_reads == 0;  # no indel
        my $indel_frac_pos = ($in_reads + $del_reads) / $n_bases;

        my @out;
        push @out, ($ref, $pos, $refcall, $cvg, "indel-target");
        push @out, ($indel_frac_pos >= $o_indelfrac and 
                    scalar(keys %indel_lengths) == 1 and
                    scalar(keys %indel_ops) == 1) ? "good" : "bad";
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
                             '@0 indel-targets -f FILE ...'
    -t | --targets TARGETS   List of indel targets to apply, output from
                             '@0 indel-targets -f FILE ...'

    --help, -?               help message

";

    # options for 

    my $o_fasta;
    my $o_targets;
    my $o_stdin;

    GetOptions(""          => \$o_stdin,
               "fasta=s"   => \$o_fasta,
               "targets=s" => \$o_targets
    ) or print_usage_and_exit($usage, "unknown option");

    open(TARGETS, "<$o_targets") or die "Could not open targets file '$o_targets': $!";
    my $target_assembly = <TARGETS>;  # first line is '#assembly:filename'
    die "First line of '$o_targets' does not contain assembly filename" if $target_ref !~ m/^#assembly:/;
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

