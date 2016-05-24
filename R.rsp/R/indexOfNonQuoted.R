###########################################################################/**
# @RdocDefault indexOfNonQuoted
#
# @title "Gets the first index of a string that is not inside a double qouted string"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{str}{The @character string to be scanned.}
#   \item{pattern}{The @character string to be searched for.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer giving the position of (the first character of)
#   the search string in the main string.  If not found, -1 is returned.
# }
#
# @author
#
# \seealso{
#   @see "base::grep".
# }
#
# @keyword programming
# @keyword utilities
# @keyword internal
#*/###########################################################################
setMethodS3("indexOfNonQuoted", "default", function(str, pattern, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'str':
  str <- as.character(str);

  # Argument 'pattern':
  pattern <- as.character(pattern);


  totalPos <- 0L;     # The position from the start of the string
  len <- 0L;          # The default match length
  qm <- NULL;         # The current qoutation mark of a string, if exists.
  ready <- FALSE;
  while(!ready) {
    # Get the first occurance of pattern in buffer
    pos <- regexpr(pattern, str);
    if (pos == -1L) return(-1L);

    totalPos <- totalPos + pos;

    tmp <- substring(str, first=1L, last=pos-1L);

    # a. Remove all espaced (doubled) backslashes, i.e. '\\'
    #    (A backslash has to be escaped in C gsub(), i.e. '\\'. Each of
    #    these two backslashes has in turn to be escaped in R doubleing
    #    the number of backslashes again!)
    tmp <- gsub("\\\\\\\\", "", tmp);

    # b. Remove all espaced quotes, i.e. '\"'.
    tmp <- gsub("\\\\[\"\']", "", tmp);

    # c. Remove all non quotes
    tmp <- gsub("[^'\"]", "", tmp);

    # d. Exclude all single or double quoted strings.
    while (nchar(tmp) > 0L) {
      if (is.null(qm)) {
        # d. Get first quotation mark
        qm <- substring(tmp, first=1L, last=1L);
        tmp <- substring(tmp, first=2L);
      }

      # e. Exclude first (single or double) quoted string.
      if (qm == "'") {
        pattern <- "^[^']*'";
      } else {
        pattern <- "^[^\"]*\"";
      }

      if (regexpr(pattern, str) != -1L) {
        str <- gsub(pattern, "", str);
        qm <- NULL;
      }
    } # while (...)

    len <- attr(pos, "match.length");
    str <- substring(str, first=pos+len);

    ready <- is.null(qm);
  } # while (!ready)

  # The found position
  pos <- as.integer(totalPos);
  attr(pos, "match.length") <- len;

  pos;
}, protected=TRUE) # indexOfNonQuoted()


##############################################################################
# HISTORY:
# 2007-04-07
# o Replaced gsub pattern "\\\[\"\']" with "\\\\[\"\']" in indexOfNonQuoted().
# 2006-07-04
# o Made the Rdoc help internal.
# 2005-09-16
# o Made the method protected.
# 2005-08-14
# o BUG FIX: Forgot to deal with single quotation marks.
# o BUG FIX: Internal countQoutationMarks() would incorrectly consider '\\"'
#   as an escaped quotation mark. Fix was to remove all '\\' first.
# 2005-08-13
# o Created. Made a seperate function because it most likely will be useful
#   elsewhere.
##############################################################################
