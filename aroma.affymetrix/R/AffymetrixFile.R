###########################################################################/**
# @RdocClass AffymetrixFile
#
# @title "The abstract AffymetrixFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFile object represents a single Affymetrix file,
#  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the generic
#  \code{fromFile()} method.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of an @see "AffymetrixFileSet".
# }
#*/###########################################################################
setConstructorS3("AffymetrixFile", function(...) {
  extend(AromaMicroarrayDataFile(...), c("AffymetrixFile",
                                             uses("AromaPlatformInterface")));
}, abstract=TRUE)



############################################################################
# HISTORY:
# 2008-05-18
# o Now "provides" AffymetrixPlatform methods.
# 2008-05-09
# o Added getPlatform().
# 2007-09-16
# o Removed obsolete copyFile() from AffymetrixFile.  Use copyTo() of
#   GenericDataFile instead.
# 2007-09-14
# O Now AffymetrixFile inherits from GenericDataFile.
# 2007-09-13
# o Added missing setAttributesByTags() to AffymetrixFile.
# 2007-08-09
# o Added static renameToUpperCaseExt().
# 2007-03-20
# o Added getAlias() and setAlias().  Note, getName() etc are still
#   unaffected by these.
# 2007-03-05
# o Added setAttributesByTags(), which now also tries to coerce values.
# o Added support for (in-memory) attributes.
# 2007-02-07
# o Added getChecksum(), writeChecksum(), readChecksum(), and
#   compareChecksum() and validateChecksum(). I did this because I noticed
#   by chance that some of my CEL files transferred via an external HDD got
#   corrupt probe signals.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o Added hasTags() and hasTag().
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-09-15
# o Added stextSize().
# 2006-08-27
# o Added stextLabel() and stextLabels(). stext is for "side text", cf.
#   mtext for "margin text". stext() is slightly more convenient than mtext
#   when it comes to different font sizes.
# o Added copyTo().
# 2006-08-14
# o Added abstract fromFile().
# 2006-08-11
# o Created from AffymetrixDataFile in order to represent CDF files too.
############################################################################
