###########################################################################/**
# @RdocDocumentation "Non-documented objects"
#
# % Other missing docs
# @eval "t <- readLines('../incl/999.missingdocs.txt'); t <- trim(unlist(strsplit(t, split=' '))); t <- t[nchar(t) > 0]; t2 <- gsub('\\[', '\\\\[', t); t <- unique(t); t <- sprintf('\\alias{%s}', t); paste(t, collapse='\n')"
#
# \description{
#   This page contains aliases for all "non-documented" objects that 
#   \code{R CMD check} detects in this package. 
#
#   Almost all of them are \emph{generic} functions that have specific 
#   document for the corresponding method coupled to a specific class. 
#   Other functions are re-defined by \code{setMethodS3()} to 
#   \emph{default} methods. Neither of these two classes are non-documented
#   in reality.
#   The rest are deprecated methods.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################

############################################################################
# HISTORY:
# 2005-05-15
# o Created to please R CMD check.
############################################################################
