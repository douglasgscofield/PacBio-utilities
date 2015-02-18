PacBio-utilities
======

A collection of scripts useful for working with PacBio-based assemblies.

The `pacbio-util` script requires [`BioPerl`][BioPerl] and [`samtools`][samtools].

[BioPerl]:  http://www.bioperl.org/
[samtools]: http://www.htslib.org/


Correction of small indels in polished assemblies
------

We have noticed that PacBio assemblies can show single-base deletions, in most
cases, in proximity to even fairly short homopolymer runs.  This does not occur
at all homopolymer runs, and in our limited samples incidence rate is perhaps 1
in 10kbp.  Single-base insertions and longer insertions may occur but are much
rarer.

There are three steps to applying this correction.  First, map a set of
Illumina reads to the PacBio assembly to generate sorted BAM file(s).  Reads
from multiple BAMs will be combined.  These reads should be from the same
individual/strain from which the assembly was created.  Reads from other
sequencing technologies, or from another individual/strain, will reduce the
effectiveness of the correction and may introduce further errors.  In all
subsequent steps it is required that the order of the sequences in the assembly
and the BAMs is identical.

Second, use the `pacbio-util` script here to generate a set of indel targets
for correction.  The command is `pacbio-util indel-targets`, and the criteria
for determining targets may be modified via options.  Use the `-h` option to
see available options and current defaults.

The third and final step is to use the command `pacbio-util indel-apply` to
apply the targets generated in the second step to the assembly to create a
corrected assembly.  Options may also be used here to restrict the set of
targets applied; use the `-h` option to see available options and current
defaults.

For an assembly `pacbio_assembly.fasta`, and mappings of Illumina paired-end
and single-end reads to this assembly in BAM files `pe.sorted.bam` and
`se.sorted.bam`, the following will generated the list of indel targets to
`targets.txt` and a corrected assembly to `corrected_assembly.fasta`.

```bash
pacbio-util indel-targets -f pacbio_assembly.fasta pe.sorted.bam se.sorted.bam > targets.txt
pacbio-util indel-apply -f pacbio_assembly.fasta -t targets.txt > corrected_assembly.fasta
```

Both commands can be combined into a one-line piped command:

```bash
pacbio-util indel-targets -f pacbio_assembly.fasta pe.sorted.bam se.sorted.bam \
| tee targets.txt \
| pacbio-util indel-apply > corrected_assembly.fasta
```

Here is example tool output (with `-v`) and scaffold size results before and
after correction.

~~~~
$ pacbio-util indel-targets -v -f pacbio_assembly.fasta pe.sorted.bam se.sorted.bam | tee targets.txt | pacbio-util indel-apply > corrected_assembly.fasta
pipe: samtools mpileup -s -B -d 10000 -L 10000 --ff 1280 -q 1 -Q 13 -f pacbio_assembly.fasta pe.sorted.bam se.sorted.bam |
[mpileup] 2 samples in 2 input files
Merging columns from 2 BAM files ...
pacbio-util indel-targets: 1281 targets generated for assembly pacbio_assembly.fasta
~~~~

Scaffold | bp before | bp after
---- | ---- | ----
unitig_1 | 6250549 | 6250752
unitig_10 | 1811588 | 1811655
unitig_7 | 15299 | 15299
unitig_52 | 7932 | 7932
unitig_46 | 8606 | 8606
unitig_14 | 1884 | 1884
unitig_0 | 10120325 | 10120600
unitig_8 | 3934972 | 3935082
unitig_12 | 240123 | 240143
unitig_3 | 16252 | 16253
unitig_54 | 22065 | 22065
unitig_44 | 12520 | 12520
unitig_53 | 5951772 | 5951911
unitig_13 | 86224 | 86224
unitig_56 | 30475 | 30475
unitig_55 | 15275 | 15275
unitig_41 | 19231 | 19232
unitig_43 | 9344 | 9344
unitig_58 | 6608738 | 6608943
unitig_9 | 2363594 | 2363680
unitig_48 | 30325 | 30325
unitig_45 | 14763 | 14763
unitig_57 | 18691 | 18685
unitig_15 | 2004 | 2004
unitig_47 | 15558 | 15558
unitig_6 | 4486498 | 4486633
unitig_60 | 26006 | 26006
unitig_61 | 16086 | 16086
unitig_51 | 9410 | 9410
unitig_59 | 14740 | 14740
unitig_4 | 19898 | 19898


### `pacbio-util indel-targets`

Generate indel targets for correction.  Uses [`samtools mpileup`][samtools] for
indel detection.  Targets are written to standard output.

Targets will only be generated where all read mappings containing an indel
agree on the size and sequence of the indel.  Further criteria for minimum
indel fraction, minimum coverage, and maximum indel size may be adjusted via
options.

~~~~
pacbio-util indel-targets -f FASTA FILE1.bam [ FILE2.bam ... ]

OPTIONS

    -f FILE, --fasta FILE    PacBio assembly in Fasta format

    FILE1.bam [ FILE2.bam ]  BAM files of reads mapped to

    -                        Read mpileup from stdin rather than running samtools.
                             Note that mapping qualities (option `-s`) are expected
                             in the mpileup input.

    --indel-frac FLOAT       Do not report indels at a position if the fraction of
                             reads containing them is below FLOAT [default 0.8]
    --indel-min-cov INT      Minimum read coverage to consider an indel [default 10]
    --include-bad-indels     Include indels in target set that do not pass our criteria,
                             these lines are marked 'bad' in the sixth column
    --max-indel-size INT     Maximum indel size to accept, beyond this the position is
                             considered in error [default 1000]
    --all-reads              Include all reads.  Default is to exclude multiply-mapped
                             reads (mapQ < 1), alignments that are not the primary
                             (SAM flag 256), and duplicate alignments (flag 1024).
    --all-bases              Include all bases.  Default is to exclude bases below
                             quality 13.
    --samtools FILE          Location of samtools executable [default samtools]

    --help, -?               help message

~~~~

Below is an example subset of targets generated with this command.  The first
line of the targets output is the name of the assembly Fasta file against which
the targets were generated.  Following that, the columns are scaffold name,
target position, reference base, total read coverage, type of target (always
`indel` for this command), whether the target passed quality criteria (`good`
or `bad`), the frequency of coverage supporting the indel, the size of the
indel (positive for insertion, negative for deletion), a string describing the
indel size and sequence, and the number of reads supporting the indel.

~~~~
#assembly:pacbio_assembly.fasta
unitig_1	101522	C	72	indel	good	0.9444	1	+1A	68
unitig_1	119319	G	127	indel	good	0.8110	1	+1A	103
unitig_1	120891	A	104	indel	good	0.9231	1	+1T	96
unitig_1	122801	A	93	indel	good	0.9032	1	+1C	84
unitig_1	125303	G	84	indel	good	0.9286	1	+1A	78
unitig_1	131734	T	67	indel	good	0.8657	1	+1A	58
unitig_1	136434	A	28	indel	good	0.8214	1	+1C	23
unitig_1	170949	A	40	indel	good	0.9500	1	+1G	38
unitig_1	171182	A	80	indel	good	0.9500	1	+1G	76
unitig_1	172016	A	67	indel	good	0.9403	1	+1C	63
unitig_1	190068	A	18	indel	good	1.0000	1	+1G	18
unitig_1	190154	T	11	indel	good	0.9091	1	+1G	10
unitig_1	191027	C	104	indel	good	0.9231	1	+1G	96
unitig_1	191173	C	96	indel	good	0.9167	1	+1T	88
unitig_1	230981	T	63	indel	good	0.9683	1	+1C	61
unitig_1	270949	G	100	indel	good	0.9800	1	+1C	98
unitig_1	307836	A	38	indel	good	0.8947	1	+1C	34
unitig_1	310152	C	105	indel	good	0.8857	1	+1A	93
unitig_1	326536	G	105	indel	good	0.9143	1	+1A	96
unitig_1	343612	G	10	indel	good	0.9000	-6	-6TTTTGT	9
unitig_1	360786	G	39	indel	good	0.8205	1	+1A	32
unitig_1	369367	C	98	indel	good	0.9592	1	+1A	94
unitig_1	371657	T	48	indel	good	0.9375	1	+1G	45
unitig_1	371776	G	43	indel	good	1.0000	1	+1T	43
unitig_1	373159	C	113	indel	good	0.9469	1	+1A	107
~~~~


### `pacbio-util indel-apply`

Apply indel targets to correct an assembly.  The corrected assembly is written
to standard output.  The order of sequences in the assembly and targets must be
the same.


~~~~
pacbio-util indel-apply [ -f FASTA ] [ -t TARGETS ]

OPTIONS

    -f | --fasta FILE        PacBio assembly in Fasta format, to be corrected.
                             Must be the same assembly used for
                             './pacbio-util indel-targets -f FILE ...'.
                             If the assembly is not specified, it is determined
                             from the first line of the target list.
    -t | --targets TARGETS   List of indel targets to apply, output from
                             './pacbio-util indel-targets -f FILE ...'
                             If no targets are specified, they are read from stdin.

    --include-bad-indels     Apply indels in target set that are marked 'bad'
    --max-indel-size INT     Maximum indel size to apply [default 1000]
    --min-indel-frac FLOAT   Minimim indel fraction to apply, a value of 0 means
                             apply all indels in target set [default 0]
    --skip-deletion-verify   Do not verify that a deletion matches the sequence at
                             that position in the assembly [default 0]
    --skip-softmasked        Do not apply targets to softmasked regions of the
                             assembly, as determined by use of lowercase [default 0]

    --help, -?               help message
~~~~

