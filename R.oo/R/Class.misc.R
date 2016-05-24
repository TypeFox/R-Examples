###########################################################################/**
# @set "class=Class"
#
# @RdocMethod getRdDeclaration
#
# @title "Gets the class declaraction in Rd format"
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
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdDeclaration", "Class", function(this, ...) {
  s <- "public"; # visibility(this);
  if (isAbstract(this))
    s <- paste(s, "abstract");
  if (isStatic(this))
    s <- paste(s, "static");
  if (inherits(this, "Class"))
    s <- paste(s, "class")
  else
    throw(getName(this), " is neither a class nor an interface.");

  s <- paste(s, " \\bold{", getName(this), "}\\cr\n", sep="");
  links <- getSuperclasses(this);

  if (length(links) > 0) {
    name <- links[1];
    link <- name;
    cls <- .getClassByName(name, mustExist=FALSE);
    if (inherits(cls, "Class")) {
      pkg <- getPackage(cls);
      if (is.null(pkg))
        link <- paste("\\link{", link ,"}", sep="")
      else
        link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
      if (isAbstract(cls))
        link <- paste("\\emph{", link, "}", sep="");
    }
    paste("\\code{", link ,"}", sep="");
    s <- paste(s, "extends ", link, "\\cr\n", sep="");
  }
  s;
}, private=TRUE);




###########################################################################/**
# @RdocMethod getRdMethods
#
# @title "Gets the methods of a class in Rd format"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{visibilities}{A @character string specifying what types of methods
#     to return.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdMethods", "Class", function(class, visibilities=c("private", "protected", "public"), ...) {
  className <- getName(class);
  methods <- getMethods(class, private=TRUE);  # Excludes empty classes!
  methods <- methods[[className]];
  methods <- names(methods);
  src <- "\\bold{Methods:}\\cr\n";

  tmpsrc <- "\\tabular{rll}{\n";
  count <- 0;
  for (method in methods) {
    fcnName <- paste(method, className, sep=".");
    fcn <- .getS3Method(fcnName);
    modifiers <- attr(fcn, "modifiers");
    if (Rdoc$isVisible(modifiers, visibilities)) {
      helpName <- Rdoc$createName(getName(class), method, escape=TRUE);
      label <- method;
      title <- Rdoc$getRdTitle(class, method);
      package <- attr(title, "package");
      if (is.null(package))
        package <- Rdoc$package;

      # Is there a specific help document for this method or not?
      if (!is.null(title)) {
  	link <- paste("\\link[", package, ":", helpName, "]{", label, "}", sep="");
      } else {
  	link <- label;
      }
      item <- paste(" \\tab \\code{", link, "} \\tab ", sep="");

      # Create the title
      if (!is.null(title)) {
  	if (title != "")
  	  item <- paste(item, title, ".\\cr", sep="");
      } else {
  	item <- paste(item, " -\\cr", sep="");
      }

      tmpsrc <- paste(tmpsrc, item, "\n", sep="");
      count <- count + 1;
    } # if(isVisible(...))
  }
  tmpsrc <- paste(tmpsrc, "}\n", sep=""); # end of \tabular{rll}

  if (count == 0)
    src <- paste(src, "\\emph{No methods defined}.\n", sep="")
  else
    src <- paste(src, tmpsrc, sep="");

  src;
}, private=TRUE);



###########################################################################/**
# @RdocMethod getRdHierarchy
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
# \seealso{
#   @seeclass
# }
#
# @keyword documentation
#*/###########################################################################
setMethodS3("getRdHierarchy", "Class", function(this, ...) {
  package <- getPackage(this);
  s <- paste("Package: ", package, "\\cr\n");
  what <- if (inherits(this, "Class")) "Class" else "Interface";
  s <- paste(s, "\\bold{", what, " ", getName(this), "}\\cr\n\n", sep="");
  indent <- "";
  for (extend in rev(getSuperclasses(this))) {
    link <- sapply(extend, FUN=function(name) {
#      isAbstract <- FALSE;
      link <- name;
      cls <- .getClassByName(name, mustExist=FALSE);
      if (inherits(cls, "Class")) {
        pkg <- getPackage(cls);
        if (is.null(pkg)) {
          link <- paste("\\link{", link ,"}", sep="")
        } else {
          link <- paste("\\link[", pkg, "]{", link ,"}", sep="");
        }
#       if (isAbstract(cls)) {
#         link <- paste("\\emph{", link, "}", sep="");
#         isAbstract <- TRUE;
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
  link <- paste("\\code{", getName(this), "}", sep="");
  if (isAbstract(this))
    link <- paste("\\emph{", link, "}", sep="");
  s <- paste(s, "\\code{", indent, "+--}", link, "\\cr\n\n", sep="");
  s;
}, private=TRUE);


#########################################################################
# HISTORY:
# 2014-03-30
# o BUG FIX: Now getRdDeclaration(), getRdHierarchy() and getRdMethods()
#   for Class handles also non-exported methods and Classes.
# 2006-05-29
# o Added support for visibility of getRdMethods().
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-10-22
# o BUG FIX: getRdMethods() returned empty \tabular{rll}{} if no methods
#   exist, but this gives an error in R CMD Rdconv.
# 2004-10-18
# o BUG FIX: Invalid filenames and link names are now escaped.
# 2004-10-17
# o Added Rdoc comments.
#########################################################################
