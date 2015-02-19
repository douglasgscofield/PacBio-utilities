TODO
====

o Add naive and model modes.  Naive mode does simple filtering, like we already do.  Model mode applies a probabilistic model.
o Analyse the targets files for the permissive and default runs
o Should I include the pieces of BioPerl I use? How do I get a list of the exact modules I use, to help with the above?
o Combine indel-targets and indel-apply?
x Decided (for now) not to use Bio::DB::Sam or another interface.
x We use BioPerl
x Consider embedding assembly filename in the targets file as the first line, behind #
x Convert pileVar.pl to pacbio-util
x It makes more sense to call samtools from within the script to verify that we are generating the pileup that we want, and to make the interface easier for users to manage.
x Come up with some options and common-sense defaults for when to generate targets
