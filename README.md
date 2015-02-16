PacBio-utilities
================

A collection of scripts useful for working with PacBio-based assemblies.  So
far this collection means one script.

**NOTE: This project has just started and does not yet do what it claims!**

pacbio-util.pl
--------------

Correct indels that can occur in PacBio assemblies in proximity to homopolymer
runs.

First, map a set of Illumina reads to the PacBio assembly.  These reads should
be from the same individual/strain from which the assembly was created.  Reads
from other sequencing technologies, or from another individual/strain, will
reduce the effectiveness of the correction and may introduce further errors.

`indel-targets` : generate indel targets for correction and summaries of their occurrence.

`indel-apply` : apply targets to correct PacBio assembly.

```bash
# map reads to reads.bam, then
samtools mpileup ... reads.bam | pacbio-util.pl indel-targets -f assembly.fasta > targets.txt
pacbio-util.pl indel-apply -f assembly.fasta -t targets.txt > corrected-assembly.fasta
```
