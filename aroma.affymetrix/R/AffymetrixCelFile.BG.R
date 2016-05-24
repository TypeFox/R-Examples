

############################################################################
# HISTORY:
# 2012-11-20
# o CLEANUP: Removed obsolete code from internal bgAdjustOptical() that
#   was never reached and that loaded affinities via obsolete APD files.
# 2011-03-26
# o ROBUSTNESS: Now internal bgAdjustGcrma(..., type="affinities") for
#   AffymetrixCelFile gives a more informative error message when there
#   are too few negative controls.
# 2010-10-02
# o We now use nomm=TRUE for all cases for type="affinities".  See code
#   for explanation.  This also solves the problems of using for instance
#   chip type MoEx-1_0-st-v1.
# 2010-09-29
# o ROBUSTNESS: Now bgAdjustGcrma(..., affinities=NULL) is deprecated and
#   throws an exception.
# o CLEANUP: Cleaned up bgAdjustGcrma().
# o Added more verbose output to bgAdjustGcrma() for AffymetrixCelFile.
# 2009-04-09 [MR]
# o BUG FIX: fixed discrepancy b/w aroma.affymetrix's gene specific binding
# adjustment and gcrma's GSB.adj
# 2009-03-29 [MR]
# o Made slight modifications for bgAdjustGcRma() to work with the
#   newer Gene 1.0 ST arrays.
# 2007-12-08
# o Now bgAdjustRma() no longer assumes 'affy' is loaded.
# 2007-06-30
# o Added .deprecated=TRUE to all methods.
# 2007-03-26
# o Added verbose output about estimated background.
# 2007-03-23
# o Replaced all usage of copyCel() with createFrom().  This allow us to
#   later update createFrom() to create CEL files of different versions.
# 2007-03-22
# o rename gsbParameters to parametersGsb to avoid clash of arguments
#   in bgAdjustGcrma.AffymetrixCelFile().  Not sure why gsbAdjust and
#   gsbParameters are being matched, but there you go.
# 2006-10-10
# o add RMA background correction
# 2006-10-06
# o make sure cdf association is inherited
# 2006-10-04
# o Debugged, tested for consistency with bg.adjust.gcrma(), docs added
# 2006-09-28
# o Created (based on AffymetrixCelFile.normalizeQuantile.R)
############################################################################
