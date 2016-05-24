###########################################################################/**
# @RdocClass DevEvalProduct
#
# @title "The DevEvalProduct class"
#
# \description{
#  @classhierarchy
#
#  A DevEvalProduct represents a handle to the "product" returned by
#  @see "devEval".
# }
#
# @synopsis
#
# \arguments{
#   \item{name, tags}{The name and optional tags of the product.}
#   \item{type}{The device type.}
#   \item{...}{Not used.}
# }
#
# \section{Fields}{
#  The following (virtual; calculate on-the-fly) fields are available:
#  \itemize{
#   \item \code{fullname}: the fullname of an image, e.g. 'foo,a,b'
#   \item \code{name}: the part of the fullname before the first comma, e.g. 'foo'
#   \item \code{tags}: the part of the fullname after the first comma, e.g. 'a,b'
#  }
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("DevEvalProduct", function(name=NULL, tags=NULL, type=NULL, ...) {
  if (!is.null(name)) {
    name <- Arguments$getCharacter(name);
    tags <- Arguments$getCharacters(tags);
    fullname <- paste(c(name, tags), collapse=",");
  } else {
    fullname <- NA_character_;
  }

  extend(BasicObject(fullname), "DevEvalProduct",
    type = type
  );
})

###########################################################################/**
# @RdocMethod as.character
# @alias as.character.DevEvalFileProduct
# @alias as.character
#
# @title "Gets a character representation of the product"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.character", "DevEvalProduct", function(x, ...) {
  getFullname(x, ...);
}, private=TRUE)


setMethodS3("view", "DevEvalProduct", function(object, ...) {})

setMethodS3("!", "DevEvalProduct", function(x) {
  view(x)
}, appendVarArgs=FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: PATCH until `[[.BasicObject` finds methods in namespaces
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("[[", "DevEvalProduct", function(x, name, ...) {
  if (name == "*") return(x);
  fcn <- sprintf("get%s", capitalize(name));
  expr <- substitute(fcn(x), list(fcn=as.name(fcn)));
  res <- tryCatch({
    eval(expr, ...);
  }, error = function(ex) {
    ex;
  });
  if (inherits(res, "Exception")) throw(res);
  if (!inherits(res, "simpleError")) return(res);
  NextMethod("[[");
}, private=TRUE)


setMethodS3("$", "DevEvalProduct", function(x, name) {
  name <- as.character(name);
  x[[name]];
}, private=TRUE)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: PATCH until `[[.BasicObject` finds methods in namespaces
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


###########################################################################/**
# @RdocMethod getFullname
# @alias getFullname.DevEvalFileProduct
# @aliasmethod getName
# @aliasmethod getTags
# @alias getFullname
# @alias getName
# @alias getTags
#
# @title "Gets the full name, name and tags"
#
# \description{
#   @get "title" consisting of a name and tags.
# }
#
# \usage{
# @usage "getFullname,DevEvalProduct"
# @usage "getFullname,DevEvalFileProduct"
# @usage "getName,DevEvalProduct"
# @usage "getTags,DevEvalProduct"
# }
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character or a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullname", "DevEvalProduct", function(this, ...) {
  as.character(unclass(this), ...);
})


setMethodS3("getName", "DevEvalProduct", function(this, ...) {
  fullname <- getFullname(this, ...);
  parts <- unlist(strsplit(fullname, split=","), use.names=FALSE);
  name <- parts[1L];
  name;
})

setMethodS3("getTags", "DevEvalProduct", function(this, collapse=",", ...) {
  fullname <- getFullname(this, ...);
  parts <- unlist(strsplit(fullname, split=","), use.names=FALSE);
  tags <- parts[-1L];
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  tags;
})

#########################################################################/**
# @RdocMethod getType
# @alias getType
#
# @title "Gets the type"
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
#  Returns a @character string or @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getType", "DevEvalProduct", function(this, ...) {
  attr(this, "type");
})




###########################################################################/**
# @RdocClass DevEvalFileProduct
#
# @title "The DevEvalFileProduct class"
#
# \description{
#  @classhierarchy
#
#  A DevEvalFileProduct is a @see "DevEvalProduct" refering to a image
#  file created by @see "devEval".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and the optional path of the image
#    file product.}
#   \item{...}{Additional arguments passed to @see "DevEvalProduct".}
# }
#
# \section{Fields}{
#  The following (virtual; calculate on-the-fly) fields are available:
#  \itemize{
#   \item \code{pathname}: the (relative) pathname of the image file, e.g. 'figures/foo,a,b.png'.  This can be used to \emph{link to} image files in for instance HTML and Markdown documents.
#   \item \code{path}, the (relative) path to the image file, e.g. 'figures/'
#   \item \code{filename}: the filename ("basename") of the image file, e.g. 'foo,a,b.png'
#   \item \code{fullname}: the fullname (the filename without the filename extension), e.g. 'foo,a,b'.  It is recommended to use this when including images in LaTeX documents that use the \code{\\usepackage{graphicx}} package.
#   \item \code{name}: the part of the fullname before the first comma, e.g. 'foo'
#   \item \code{tags}: the part of the fullname after the first comma, e.g. 'a,b'
#   \item \code{dataURI}: the Base64-encoded Data URI representation of the image, e.g. 'data:image/png;base64,iVBORw0KGgoAAAA...'.  This can be used to \emph{include} ("inlining") image files in for instance self-contained HTML and Markdown documents.
#   \item \code{data}: the character content of the image file.  This can be used to \emph{include} ("inlining") WebGL HTML-based image files in for instance self-contained HTML and Markdown documents.
#  }
# }
#
# \seealso{
#   In order to retrieve the Data URI, the \pkg{base64enc} package
#   must be installed.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("DevEvalFileProduct", function(filename=NULL, path=NULL, ...) {
  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
    path <- dirname(pathname);
    filename <- basename(pathname);
  } else {
    pathname <- NA_character_;
  }

  this <- extend(DevEvalProduct(pathname, ...), "DevEvalFileProduct");

  # Infer 'type' from pathname?
  type <- getType(this);
  if (is.null(type) && !is.na(pathname)) {
    ext <- getExtension(this);
    type <- .devTypeName(ext);
    attr(this, "type") <- type;
  }

  this;
})

setMethodS3("as.character", "DevEvalFileProduct", function(x, ...) {
  getPathname(x, ...);
}, private=TRUE)


setMethodS3("getFullname", "DevEvalFileProduct", function(this, ...) {
  filename <- getFilename(this, ...);
  gsub("[.]([^.]*)$", "", filename);
})


###########################################################################/**
# @RdocMethod getPathname
# @aliasmethod getFilename
# @aliasmethod getPath
# @aliasmethod getExtension
# @alias getPathname
# @alias getFilename
# @alias getPath
# @alias getExtension
#
# @title "Gets the (relative) pathname, filename and path"
#
# \description{
#   @get "title".
# }
#
# \usage{
# @usage "getPathname,DevEvalFileProduct"
# @usage "getFilename,DevEvalFileProduct"
# @usage "getPath,DevEvalFileProduct"
# @usage "getExtension,DevEvalFileProduct"
# }
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPathname", "DevEvalFileProduct", function(this, relative=TRUE, ...) {
  pathname <- as.character(unclass(this), ...);
  if (relative) {
    pathname <- getRelativePath(pathname);
  } else {
    pathname <- getAbsolutePath(pathname);
  }
  pathname;
})

setMethodS3("getFilename", "DevEvalFileProduct", function(this, ...) {
  basename(getPathname(this, ...));
})

setMethodS3("getPath", "DevEvalFileProduct", function(this, ...) {
  dirname(getPathname(this, ...));
})

setMethodS3("getExtension", "DevEvalFileProduct", function(this, ...) {
  filename <- getFilename(this, ...);
  gsub(".*[.]([^.]*)$", "\\1", filename);
}, private=TRUE)

setMethodS3("view", "DevEvalFileProduct", function(object, ...) {
  pathname <- object

  # WORKAROUND: browseURL('foo/bar.html', browser=NULL), which in turn
  # calls shell.exec('foo/bar.html'), does not work on Windows, because
  # the OS expects backslashes.  [Should shell.exec() convert to
  # backslashes?]  By temporarily setting the working directory to that
  # of the file, this works around this issue.
  # Borrowed from R.rsp. /HB 2014-09-19
  if (isFile(pathname)) {
    path <- dirname(pathname);
    pathname <- basename(pathname);
    opwd <- getwd();
    on.exit(setwd(opwd));
    setwd(path);
  }
  browseURL(pathname, ...)

  invisible(object)
})



#########################################################################/**
# @RdocMethod getMimeType
# @aliasmethod getMime
# @alias getMimeType
# @alias getMime
#
# @title "Gets the MIME type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{default}{The value returned, if the MIME type could not be inferred.}
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
#*/#########################################################################
setMethodS3("getMimeType", "DevEvalFileProduct", function(this, default="", ...) {
  ext <- getExtension(this, ...);
  ext <- .devTypeName(ext);
  mimeTypes <- c(
    bmp="image/bmp",
    gif="image/gif",
    jpg="image/jpeg",
    pdf="application/pdf",
    eps="application/postscript",
    postscript="application/postscript",
    png="image/png",
    svg="image/svg+xml",
    tiff="image/tiff"
  )
  mime <- mimeTypes[ext];
  if (is.na(mime)) mime <- default;
  mime;
})

setMethodS3("getMime", "DevEvalFileProduct", function(this, ...) {
  getMimeType(this, ...);
}, private=TRUE)



#########################################################################/**
# @RdocMethod getDataURI
# @alias getDataURI
# @aliasmethod getData
# @alias getData
#
# @title "Gets content as a Base64-encoded data URI"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mime}{The MIME type to be embedded in the data URI.}
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
#*/#########################################################################
setMethodS3("getDataURI", "DevEvalFileProduct", function(this, mime=getMimeType(this), ...) {
  dataURI(file=getPathname(this), mime=mime, encoding="base64");
})

setMethodS3("getData", "DevEvalFileProduct", function(this, mode=c("character", "raw"), ...) {
  # Argument 'mode':
  mode <- match.arg(mode);

  pathname <- this;
  size <- file.info(pathname)$size;
  if (mode == "character") {
    res <- readChar(con=pathname, nchars=size);
  } else if (mode == "raw") {
    res <- readBin(con=pathname, what=raw(), n=size);
  }
  res
}) # getData()


############################################################################
# HISTORY:
# 2014-09-17
# o Added view() and !() to DevEvalProduct.
# 2014-09-15
# o BUG FIX: Now getPathname(..., relative=FALSE) returns the absolute
#   pathname.
# 2014-09-02
# o Added getData() to DevEvalFileProduct.
# 2013-09-17
# o ROBUSTNESS: Now getDataURI() throws an Exception is suggested
#   package 'base64enc' is not installed.
# 2013-08-29
# o Now getPathname() returns the relative pathname, by default.
# 2013-08-27
# o Added the DevEvalProduct and DevEvalFileProduct classes.
# o Created.
############################################################################
