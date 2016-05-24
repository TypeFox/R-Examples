###########################################################################/**
# @RdocDefault rargs
#
# @title "Gets RSP arguments of an RSP document"
#
# \description{
#  @get "title", if any.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "rcompile".}
# }
#
# \value{
#   Returns a @data.frame of class 'RspArguments'.
# }
#
# \details{
#   Any RSP preprocessing variable with an 'description' attribute
#   is considered to be an RSP argument.
# }
#
# @examples "../incl/rargs.Rex"
#
# @author
#
# \seealso{
#  @see "rfile".
# }
#
# @keyword file
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("rargs", "default", function(...) {
  # Parse RSP document
  doc <- rcompile(..., until="directives", output=RspDocument());

  # Extract RSP preprocessing directives
  keep <- unlist(sapply(doc, FUN=inherits, "RspUnparsedDirective"), use.names=FALSE);
  doc <- doc[keep];

  # Parse RSP directives
  for (idx in seq_along(doc)) {
    doc[[idx]] <- parse(doc[[idx]]);
  }

  # Extract RSP preprocessing variables
  keep <- unlist(sapply(doc, FUN=inherits, "RspVariableDirective"), use.names=FALSE);
  doc <- doc[keep];

  # Subset by those with 'description' attributes.
  keep <- unlist(sapply(doc, FUN=hasAttribute, "description"), use.names=FALSE);
  doc <- doc[keep];

  args <- lapply(doc, FUN=function(d) {
    attrs <- getNameContentDefaultAttributes(d);
    type <- as.character(d);
    default <- attrs$default;
    if (is.null(default)) default <- vector(mode=type, length=1L);
    default <- as(default, type);
    data.frame(
      name        = attrs$name,
      type        = type,
      default     = default,
      description = attr(d, "description"),
      stringsAsFactors=FALSE
    )
  });
  if (length(args) > 0L) {
    args <- Reduce(rbind, args);
  } else {
    # Default
    args <- data.frame(
      name        = "",
      type        = "",
      default     = "",
      description = "",
      stringsAsFactors=FALSE
    );
    args <- args[c(),];
  }

  rownames(args) <- NULL;
  class(args) <- c("RspArguments", class(args));
  args;
}) # rargs()


#########################################################################/**
# @set "class=RspArguments"
# @RdocMethod print
#
# @title "Prints RSP arguments"
#
# \description{
#  @get "title" returned by @see "rargs".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @see "rargs".
# }
#*/#########################################################################
setMethodS3("print", "RspArguments", function(x, ...) {
  s <- NULL;
  for (kk in seq_len(nrow(x))) {
    arg <- x[kk,];
    title <- sprintf("'%s' [%s]", arg$name, arg$type);
    title <- sprintf("%s:", title);
    if (!is.null(arg$default)) {
      default <- sprintf("    Default: '%s'", arg$default);
    } else {
      default <- "    Default:";
    }
    desc <- arg$description;
    desc <- sprintf("    %s", desc);
    s <- c(s, title, default, desc, "");
  } # for (kk ...)
  s <- paste(s, collapse="\n");
  cat(s, "\n", sep="");
})


############################################################################
# HISTORY:
# 2013-04-02
# o Created.
############################################################################
