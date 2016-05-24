###########################################################################/**
# @RdocClass GenericReporter
#
# @title "The GenericReporter class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{tags}{A @character @vector of tags to be added to the output path.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("GenericReporter", function(tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getTags(tags, collapse=NULL);
  }


  extend(Object(), "GenericReporter",
    .alias = NULL,
    .tags = tags
  )
})


setMethodS3("as.character", "GenericReporter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getAlias
#
# @title "Gets the alias of the report"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character, or @NULL if no alias is set.
# }
#
# \seealso{
#   @seemethod "setAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlias", "GenericReporter", function(this, ...) {
  this$.alias;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the alias of the report"
#
# \description{
#   @get "title".
#   If specified, the alias overrides the report name, which is used by
#   default.
# }
#
# @synopsis
#
# \arguments{
#  \item{alias}{A @character string for the new alias of the report.
#   The alias must consists of valid filename characters, and must not
#   contain commas, which are used to separate tags.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \seealso{
#   @seemethod "getAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setAlias", "GenericReporter", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    alias <- Arguments$getFilename(alias);  # Valid filename?

    # Assert that no commas are used.
    if (regexpr("[,]", alias) != -1) {
      throw("Aliases (names) must not contain commas: ", alias);
    }
  }

  this$.alias <- alias;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \details{
#  If a name alias has not been set explicitly, the input name will be used.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "GenericReporter", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    name <- getInputName(this);
  }
  name;
})




###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the reporter"
#
# \description{
#  @get "title", which are the input tags plus additional tags.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "GenericReporter", function(this, collapse=NULL, ...) {
  tags <- getInputTags(this);

  tags <- c(tags, this$.tags);

  # In case this$.tags is not already split
  tags <- Arguments$getTags(tags, collapse=NULL);

  tags <- locallyUnique(tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Keep non-empty tags
  tags <- tags[nzchar(tags)];

  tags <- locallyUnique(tags);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    if (length(tags) > 0)
      tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})

setMethodS3("getInputName", "GenericReporter", abstract=TRUE, protected=TRUE)

setMethodS3("getInputTags", "GenericReporter", abstract=TRUE, protected=TRUE)

setMethodS3("getAsteriskTags", "GenericReporter", function(this, ...) {
  "";
}, protected=TRUE)


setMethodS3("getFullName", "GenericReporter", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getReportSet", "GenericReporter", abstract=TRUE, protected=TRUE);


setMethodS3("getRootPath", "GenericReporter", function(this, ...) {
  "reports";
}, protected=TRUE)


setMethodS3("getMainPath", "GenericReporter", function(this, ...) {
  # Create the (sub-)directory tree for the report set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  name <- getName(this);

  # Tags
  tags <- getTags(this, collapse=",");
  if (length(tags) == 0 || !nzchar(tags)) {
    tags <- "raw";  # Default
  }

  # The full path
  path <- filePath(rootPath, name, tags);
  path <- Arguments$getWritablePath(path);

  path;
}, protected=TRUE)



setMethodS3("getPath", "GenericReporter", abstract=TRUE);


setMethodS3("setup", "GenericReporter", abstract=TRUE);


###########################################################################/**
# @RdocMethod process
#
# @title "Generates report"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{See subclasses.}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "GenericReporter", abstract=TRUE);



##############################################################################
# HISTORY:
# 2008-04-12
# o Created from AffymetrixCelSetReporter.R.
##############################################################################
