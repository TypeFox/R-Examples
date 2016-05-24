###########################################################################/**
# @set "class=classRepresentation"
#
# @RdocMethod getKnownSubclasses
# @keyword internal
#
# @title "Gets the known subclasses"
#
# \description{
#   @get "title".
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
# @author
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getKnownSubclasses", "classRepresentation", function(this, ...) {
  this@subclasses$signature;
})



###########################################################################/**
# @RdocMethod getSuperclasses
# @keyword internal
#
# @title "Gets the superclasses"
#
# \description{
#   @get "title".
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
# @author
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getSuperclasses", "classRepresentation", function(this, ...) {
  superClasses <- NULL;
  for (contain in attr(this, "contains")$vector) {
    superClasses <- c(superClasses, contain@superClass);
  }
  superClasses;
})



###########################################################################/**
# @RdocMethod getRdHierarchy
# @keyword internal
#
# @title "Gets the class hierarchy in Rd format"
#
# \description{
#   @get "title".
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
# @author
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdHierarchy", "classRepresentation", function(this, ...) {
  package <- "???";
  name <- this@className;
  superClasses <- getSuperclasses(this);

  s <- paste("Package: ", package, "\\cr\n");
  s <- paste(s, "\\bold{Class ", name, "}\\cr\n\n", sep="");
  indent <- "";
  for (extend in rev(superClasses)) {
    link <- sapply(extend, FUN=function(name) {
      link <- name;
      if (isClass(name)) {
        cls <- getClass(name);
        link <- paste("\\link{", link ,"}", sep="")
      }
      paste("\\code{", link ,"}", sep="");
    });
    if (indent == "") {
      s <- paste(s, link, "\\cr\n", sep="");
      indent <- "~~";
    } else {
      s <- paste(s, "\\code{", indent, "+--}", link, "\\cr\n", sep="");
      indent <- paste(indent, "~~~~~", sep="");
    }
    s <- paste(s, "\\code{", indent, "|}\\cr\n", sep="");
  }
  link <- paste("\\code{", name, "}", sep="");
  s <- paste(s, "\\code{", indent, "+--}", link, "\\cr\n\n", sep="");
  s;
}, private=TRUE);






###########################################################################/**
# @RdocMethod getRdDeclaration
# @keyword internal
#
# @title "Gets the class declaration in Rd format"
#
# \description{
#   @get "title".
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
# @author
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdDeclaration", "classRepresentation", function(this, ...) {
  name <- this$className;

  s <- "public"; # visibility(this);
  s <- paste(s, "class")

  s <- paste(s, " \\bold{", name, "}\\cr\n", sep="");
  links <- getSuperclasses(this);

  if (length(links) > 0) {
    name <- links[1];
    link <- name;
    # TO DO/FIX ME: This part only works when packages are attached.
    # /HB 2013-10-08
    if (exists(name, mode="function")) {
      cls <- get(name, mode="function");
      if (inherits(cls, "Class")) {
        pkg <- getPackage(cls);
        if (is.null(pkg))
          link <- paste("\\link{", link ,"}", sep="")
        else
          link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
        if (isAbstract(cls))
          link <- paste("\\emph{", link, "}", sep="");
      }
    }
    paste("\\code{", link ,"}", sep="");
    s <- paste(s, "extends ", link, "\\cr\n", sep="");
  }
  s;
}, private=TRUE);





###########################################################################/**
# @RdocMethod getRdMethods
# @keyword internal
#
# @title "Gets the methods in Rd format"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{visibility}{-}
#   \item{trial}{-}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdMethods", "classRepresentation", function(class, visibility=c("public", "protected", "private"), trial=FALSE, ...) {
  src <- NULL;
  src <- paste(src, "\\emph{No methods defined}.\n", sep="")
  src;
}, private=TRUE);



#########################################################################
# HISTORY:
# 2005-06-08
# o Added keyword "internal" to all methods, because of change in Rdoc.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-10-17
# o Added Rdoc comments.
#########################################################################
