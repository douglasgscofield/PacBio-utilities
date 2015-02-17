TODO
====

x Convert pileVar.pl to pacbio-util
* It makes more sense to call samtools from within the script to verify that we are generating the pileup that we want, and to make the interface easier for users to manage.
* Come up with some options and common-sense defaults for when to generate targets
* Use Bio::DB::Sam?
* Should I include the pieces of BioPerl I use?
* In Perl, how do I get a list of the exact modules I use, to help with the above?
* Combine indel-targets and indel-apply?
* We use BioPerl
x Consider embedding assembly filename in the targets file as the first line, behind #
* Automate quality base detection
