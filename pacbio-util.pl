#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use File::Which;
use Getopt::Long;
use List::Util qw/sum/;

my $o_fasta;
my $o_BAM;
my $o_stdin;

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
sub target_apply();

################
# main loop

my $command = shift @ARGV;

if ($command eq "indel-targets") {
    indel_targets();
} elsif ($command eq "target-apply") {
    target_apply();
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

    GetOptions(""          => \$o_stdin,
               "fasta=s"   => \$o_fasta
    ) or print_usage_and_exit($usage, "unknown option");
    my @BAM = grep { -f or die "cannot find BAM: $_" } @ARGV;
    my $n_BAM = scalar(@BAM);

    my $samtools = which('samtools') or die "$command: samtools not found";

    my $samtools_pipe = "$samtools mpileup -s -BQ0 -d 1000000 -L 1000000 -f $o_fasta ".join(" ", @BAM)." |";
    open PILEUP, $samtools_pipe or die "Could not initiate samtools pipe: $!";

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

sub target_apply() {
    print STDERR "you've reached target_apply\n";
}

