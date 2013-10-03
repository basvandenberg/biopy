biopy
=====

The biopy package is a Swiss Army knife for bioinformaticians, albeit still
a rather small one. It is basically a collection of code utilities I used
during the past few years and therefore it contains a lot of undocumented stuff
that is probably pretty useless, but there is also a lot of useful stuff in
there.

roc.py
------

Contains functions for calculating and plotting receiver operator
characteristic (ROC) curves. The ROC class can be used for plotting single ROC
curves and the RocCollection class can be used for plotting multiple ROC-curves
in one plot. This could be useful for comparing different ROC curves or for
showing ROC curves for all cross-validation loops together with their average
ROC-curve.

sequtil.py
----------

Contains many util functions for biological sequences (DNA and protein). It
primarily contains functions for calculating sequence-based protein features,
such as amino acid composition, codon composition, pseudo-amino acid
composition, autocorrelation and so on.

file\_io.py
-----------

File parsers, most of which will not be useful for general use. It contains
a very basic FASTA file parser that might be useful.


Dependencies
============

The sofware is developed for python2.7. The dependencies for using the software
are:

- numpy >= 1.7.1
- scipy >= 0.12.0
- matplotlib >= 1.2.2


Installation
============

For installation on a Linux system, use:

    sudo python setup.py install
