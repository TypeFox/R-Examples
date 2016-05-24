#########################################################################/**
# @RdocPackage rtf
#
# \encoding{latin1}
# 
# \description{
# A set of R functions to output Rich Text Format (RTF) files with high
# resolution tables and graphics that may be edited with a standard word
# processor such as Microsoft Word.  This is useful for compiling R
# results into a document for further editing or to be merged with another
# document.
#
# While Sweave and odfWeave provide nicely formatted outputs, the syntax
# can be difficult for beginners to master. Moreover, Sweave compilation
# to a PDF is not editable.  Other packages such as SWord and R2wd provide
# similar functionality but require the user to be running a Microsoft
# Windows OS.
# 
# This package is written in pure R and does not require leaving the R
# environment to write files. R data frames and table objects are converted to
# nicely formatted RTF tables.  One important limitation of the RTF 
# specification is that vector graphics output is limited to Windows Meta File
# (WMF) and Enhanced Meta File (EMF) formats.  Because these formats are not
# supported across platforms, this package currently only supports RTF
# embedding of PNG plots and images.  To ensure high quality reports, the
# resolution may be specified when writing the RTF output.
# }
#
# \section{Requirements}{
#   This package depends on the \pkg{R.oo} package.
# }
#
# \section{Usage}{
#   For usage details @see "RTF".
# }
#
# @author
# 
# \references{
# [1] \url{http://en.wikipedia.org/wiki/Rich_Text_Format}\cr
#
# [2] \url{http://latex2rtf.sourceforge.net/rtfspec_7.html#rtfspec_paraforprop}\cr
# }
#
# \keyword{ rtf }
#
#*/#########################################################################
