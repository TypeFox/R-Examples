# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().


############################################################################
# HISTORY:
# 2013-08-03 [HB]
# o CLEANUP: Removed appendVarArgs(writeCdf).
# 2012-08-30 [HB]
# o CLEANUP: No longer a need for appendVarArgs(write).
# o CLEANUP: No longer a need for appendVarArgs(boxplot.stats).
# o CLEANUP: No longer a need for appendVarArgs(getPackageName).
# 2011-01-09 [HB]
# o Added appendVarArgs(writeCdf).
# 2010-02-08 [HB]
# o Added appendVarArgs(boxplot.stats) so that one can pass argument
#   'show.names' to bxp() via plotRle().
# 2007-02-27 [HB]
# o BUG FIX: Removed explicit reference to 'base' etc again. The reason is
#   that if a previous package already modified, say, write(), to become a
#   generic function, that was overwritten again when this package was
#   loaded.
# 2007-02-23 [KS]
# o Make explicit reference to 'base' - this is safer, in case of colMeans()
#   defined by user or other packages.
# 2006-03-24
# o Created to please R CMD check.
############################################################################
