#########################################################################/**
# @set "class=Arguments"
# @RdocMethod getTags
#
# @title "Gets and validates tags"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{@character @vector of tags.}
#   \item{na.rm}{If @TRUE, empty ("missing") tags are dropped.}
#   \item{collapse}{A @character string specifying how the tags should
#     be concatenated into a single string.
#     If @NULL, they are not concattenated.}
# }
#
# \value{
#  Returns a @character string or
#  @character @vector (iff \code{collapse} is @NULL).
# }
#
# @author
#
# \seealso{
#   For more information see \code{\link[R.utils]{Arguments}}.
# }
#*/#########################################################################
setMethodS3("getTags", "Arguments", function(static, ..., na.rm=TRUE, collapse=",") {
  args <- list(...)

  ## Drop empty elements
  keep <- (unlist(lapply(args, FUN=length), use.names=FALSE) > 0)
  args <- args[keep]

  ## Nothing to do?
  if (length(args) == 0) return(NULL)

  ## Drop missing elements?
  if (na.rm) {
    keep <- !unlist(lapply(args, FUN=is.na), use.names=FALSE)
    args <- args[keep]
  }

  # Generate tags
  tags <- do.call(paste, args=c(args, sep=",", collapse=","))
  tags <- Arguments$getCharacters(tags);
  tags <- strsplit(tags, split=",", fixed=TRUE);
  tags <- unlist(tags);
  tags <- trim(tags);

  # Drop missing tags?
  if (na.rm) {
    tags <- tags[!is.na(tags)];
  }

  # Drop empty tags
  tags <- tags[nchar(tags, type="chars") > 0L];

  # Nothing to do?
  if (length(tags) == 0) return(NULL);

  # Collapse?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2012-10-16
# o Added Rdoc comments.
# 2011-02-18
# o Moved getTags() for Arguments from aroma.core to R.filesets.
# 2010-01-25
# o Added static getTags() to Arguments.
# ...
############################################################################
