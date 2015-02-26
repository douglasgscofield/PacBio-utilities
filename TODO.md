TODO
====

- Add filter and model modes.  Filter mode does naive filtering, like we already do.  Model mode applies a probabilistic model, and is tuned for these specific situations.
- Analyse the targets files for the permissive and default runs
- Should I include the pieces of BioPerl I use?
- How do I get a list of the exact modules I use, to help with the above?
- Combine indel-targets and indel-apply?

Completed
---------
* Decided (for now) not to use Bio::DB::Sam or another interface.
* We use BioPerl
* Consider embedding assembly filename in the targets file as the first line, behind #
* Convert pileVar.pl to pacbio-util
* It makes more sense to call samtools from within the script to verify that we are generating the pileup that we want, and to make the interface easier for users to manage.
* Come up with some options and common-sense defaults for when to generate targets
