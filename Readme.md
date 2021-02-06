# Pindel

This is a copy of an applet created by DNAnexus, which implements the Pindel tool (version 0.2.5a7) for detection of structural variants (large deletions, medium sized insertions, inversions, tandem duplications, etc).

The following is the Readme as provided by DNAnexus and links are now out of date. Use https://github.com/genome/pindel/ instead.

## What does this app do?

Pindel finds breakpoints for deletions, insertions, inversions, and
duplications. It detects these on a base pair level.

## What are typical use cases for this app?

This app can be run on whole genome or whole exome data to try and find
base pair level calls for structural variants.

## What data are required for this app to run?

This app requires an input file containing the mappingsand a reference genome
in fasta format (`*.fa`). 

This app must have one of the following input combinations or it will raise an AppError:
1.  BAM file(s) with extension '.bam' from BWA/MOSAIK AND BAM config file (as specied at http://gmt.genome.wustl.edu/pindel/0.2.4/user-manual.html) OR insert size
2.  BAM file(s) with extension '.bam' not from BWA/MOSAIK AND bam_not_produced_by_bwa=True (default=False) AND sequence_platform AND insert size
3.  Pindel input file (as specified at http://gmt.genome.wustl.edu/pindel/0.2.4/user-manual.html)


## What does this app output?

This app outputs 6 files which each contain a different piece of structural
variation information.

A _BP file which contains breakpoints not assigned to any event.
A _D file which contains deletions.
A _INV file which contains inversions.
A _LI file which contains long insertions.
A _SI file which contains short insertions.
A _TD file which contains tandem duplications

This app can also output a VCF containing some of the structural variant
information. Because the app includes the reads in its file, the files
can become very large.

## How does this app work?
This app uses a pattern growth approach to find the breaakpoints of
variants from paired-end short reads.

For more informatoin, consult the Pindel website at:

http://gmt.genome.wustl.edu/pindel/0.2.4/index.html
