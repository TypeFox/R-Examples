###########################################################################/**
# @RdocClass RspFileProduct
#
# @title "The RspFileProduct class"
#
# \description{
#  @classhierarchy
#
#  An RspFileProduct is an @see RspProduct that represents an
#  RSP product in form of a file, e.g. LaTeX, Sweave and knitr documents.
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{An existing file.}
#   \item{...}{Additional arguments passed to @see "RspProduct".}
#   \item{mustExist}{If @TRUE, it is asserted that the file exists.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspFileProduct", function(pathname=NA, ..., mustExist=TRUE) {
  # Argument 'pathname':
  if (!is.null(pathname) && !is.na(pathname)) {
    if (!isUrl(pathname)) {
      withoutGString({
        pathname <- Arguments$getReadablePathname(pathname, mustExist=mustExist);
      })
    }
  }

  extend(RspProduct(pathname, ...), "RspFileProduct");
})


setMethodS3("print", "RspFileProduct", function(x, ...) {
  s <- sprintf("%s:", class(x)[1L]);

  s <- c(s, sprintf("Pathname: %s", x));

  # File size
  fileSize <- getFileSize(x, "units");
  if (!is.na(fileSize)) {
    fileSizeB <- sprintf("%.0f bytes", getFileSize(x, "numeric"));
    if (fileSizeB != fileSize) {
      fileSize <- sprintf("%s (%s)", fileSize, fileSizeB);
    }
  }
  s <- c(s, sprintf("File size: %s", fileSize));

  s <- c(s, sprintf("Content type: %s", getType(x)));

  md <- getMetadata(x, local=FALSE);
  for (key in names(md)) {
    s <- c(s, sprintf("Metadata '%s': '%s'", key, md[[key]]));
  }

  s <- c(s, sprintf("Has processor: %s", hasProcessor(x)));

  s <- paste(s, collapse="\n");
  cat(s, "\n", sep="");
}, protected=TRUE)



setMethodS3("view", "RspFileProduct", function(object, ...) {
  # WORKAROUND: browseURL('foo/bar.html', browser=NULL), which in turn
  # calls shell.exec('foo/bar.html'), does not work on Windows, because
  # the OS expects backslashes.  [Should shell.exec() convert to
  # backslashes?]  By temporarily setting the working directory to that
  # of the file, view() for RspFileProduct works around this issue.
  if (isFile(object)) {
    path <- dirname(object);
    pathname <- basename(object);
    opwd <- getwd();
    on.exit(setwd(opwd));
    setwd(path);
  } else {
    pathname <- object;
  }
  browseURL(pathname, ...);
  invisible(object);
}, proctected=TRUE)



setMethodS3("getType", "RspFileProduct", function(object, default=NA_character_, as=c("text", "IMT"), ...) {
  as <- match.arg(as);
  res <- NextMethod("getType", default=NA_character_);

  if (is.na(res)) {
    # Infer type from the filename extension?
    if (isFile(object) || isUrl(object)) {
      res <- extensionToIMT(object);
    }
  }

  # Fall back to a default?
  if (is.na(res)) {
    default <- as.character(default);
    res <- default;
  }

  if (as == "IMT" && !is.na(res)) {
    res <- parseInternetMediaType(res);
  }

  res;
}, protected=TRUE)


setMethodS3("getFileSize", "RspFileProduct", function(object, what=c("numeric", "units"), sep="", ...) {
  # Argument 'what':
  what <- match.arg(what);

  pathname <- object
  if (is.null(pathname) && isUrl(pathname)) {
    fileSize <- NA_real_;
  } else {
    fileSize <- file.info2(pathname)$size;
  }

  if (what == "numeric")
    return(fileSize);

  if (is.na(fileSize))
    return(fileSize);

  units <- c("bytes", "kB", "MB", "GB", "TB");
  scale <- 1;
  for (kk in seq_along(units)) {
    unit <- units[kk];
    if (fileSize < 1000)
      break;
    fileSize <- fileSize/1024;
  }
  fileSize <- sprintf("%.2f %s%s", fileSize, sep, unit);
  fileSize <- gsub(".00 bytes", " bytes", fileSize, fixed=TRUE);

  fileSize;
})



setMethodS3("findProcessor", "RspFileProduct", function(object, ..., verbose=FALSE) {
  isFALSE <- function(x) {
    if (is.character(x)) return(toupper(x) == "FALSE")
    if (is.logical(x) && !x) return(TRUE)
    FALSE
  }

  localCompileLaTeX <- function(..., texinputs=NULL) {
    if (!is.null(source)) {
      path <- dirname(source)
      pathA <- getAbsolutePath(path)
      texinputs <- c(texinputs, pathA)
    }
    compileLaTeX(..., texinputs=texinputs)
  } # localCompileLaTeX()


  localCompressPDF <- function(pathname, ..., verbose=FALSE) {
    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
      pushState(verbose)
      on.exit(popState(verbose))
    }

    ## Disabling further postprocessing after this,
    ## which avoids recursive loop.
    metadata <- attr(pathname, "metadata")
    metadata$postprocess <- FALSE
    attr(pathname, "metadata") <- metadata

    ## Get compression level
    compression <- metadata$compression

    verbose && enter(verbose, "Trying to compress PDF")
    verbose && cat(verbose, "Compression: ", compression)

    pathT <- tempfile(pattern=".dir", tmpdir=".")
    on.exit(removeDirectory(pathT, recursive=TRUE), add=TRUE)

    ## ROBUSTNESS: If compression fails for one reason or the
    ## other, fall back to keep the non-compressed version.
    tryCatch({
      verbose && enter(verbose, "R.utils::compressPDF()...")
      suppressWarnings({
        pathnameZ <- compressPDF(pathname, outPath=pathT,
                                 compression=compression)
      })
      verbose && exit(verbose)

      size <- file.info(pathname)$size
      sizeZ <- file.info(pathnameZ)$size
      if (!identical(sizeZ, size)) {
        renameFile(pathnameZ, pathname, overwrite=TRUE)
      }
    }, error = function(ex) {
      msg <- sprintf("Compression of '%s' using '%s' failed. Keeping the original PDF file. Reason was: %s", pathname, compression, ex$message)
      verbose && cat(verbose, msg)
      warning(msg)
    })

    verbose && exit(verbose)

    pathname
  } # localCompressPDF()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating document-type specific processor");
  type <- getType(object);
  verbose && cat(verbose, "RSP product content type: ", type);

  # Nothing to do?
  if (is.na(type)) {
    verbose && cat(verbose, "Processor found: <none>");
    verbose && exit(verbose);
    return(NULL);
  }
  type <- parseInternetMediaType(type)$contentType;


  # Nothing to do?
  postprocess <- getMetadata(object, "postprocess", local=TRUE)
  if (isFALSE(postprocess)) {
    verbose && cat(verbose, "Processing disabled: metadata variable 'postprocess' is FALSE")
    verbose && exit(verbose)
    return(NULL)
  }


  source <- getMetadata(object, "source", local=TRUE);
  if (is.null(source)) {
    verbose && cat(verbose, "Source document: <unknown>");
  } else {
    verbose && cat(verbose, "Source document: ", sQuote(source));
  }



  # Find a down-stream compiler/processor:
  fcn <- switch(type,
    # RSP documents:
    # *<ext>.rsp => *.<ext>
    "application/x-rsp" = function(...) { compileRsp(..., postprocess=FALSE) },

    # LaTeX documents:
    # *.tex => ... => *.pdf
    "application/x-tex" = localCompileLaTeX,
    "application/x-latex" = localCompileLaTeX,

    ## PDF documents:
    # *.pdf => ... => *.pdf
    "application/pdf" = localCompressPDF,

    # Markdown documents:
    # *.md => *.html
    "application/x-markdown" = compileMarkdown,

    # Markdown documents:
    # *.txt => *.html, ...
    "application/x-asciidoc" = function(...) { compileAsciiDoc(..., postprocess=FALSE) },

    # Sweave Rnw documents:
    # *.Rnw => *.tex
    "application/x-sweave" = function(...) { compileSweave(..., postprocess=FALSE) },

    # Knitr Rnw documents:
    # *.Rnw => *.tex
    "application/x-knitr" = function(...) { compileKnitr(..., postprocess=FALSE) },

    # Knitr Rmd documents:
    # *.Rmd => *.html
    "application/x-rmd" = function(...) { compileKnitr(..., postprocess=FALSE) },
    # Knitr Rhtml documents:
    # *.Rhtml => *.html
    "application/x-rhtml" = function(...) { compileKnitr(..., postprocess=FALSE) },

    # Knitr Rtex documents:
    # *.Rtex => *.tex
    "application/x-rtex" = function(...) { compileKnitr(..., postprocess=FALSE) },

    # Knitr Rrst documents:
    # *.Rrst => *.rst
    "application/x-rrst" = function(...) { compileKnitr(..., postprocess=FALSE) },

    # AsciiDoc Rnw documents:
    # *.Rnw => *.txt
    "application/x-asciidoc-noweb" = function(...) { compileAsciiDocNoweb(..., postprocess=FALSE) },

    # Sweave or Knitr Rnw documents:
    # *.Rnw => *.tex
    "application/x-rnw" = function(...) { compileRnw(..., postprocess=FALSE) }
  );

  if (is.null(fcn)) {
    verbose && cat(verbose, "Processor found: <none>");
  } else {
    # Get the metadata attributes
    metadata <- getMetadata(object, local=TRUE);

    # Make sure the processor returns an RspFileProduct
    fcnT <- fcn
    processor <- function(...) {
       do.call(fcnT, args=c(list(...), list(metadata=metadata)));
    }

    fcn <- function(pathname, ...) {
      # Arguments 'pathname':
      if (!isUrl(pathname)) {
        withoutGString({
          pathnameT <- Arguments$getReadablePathname(pathname);
          pathname[1] <- pathnameT;  ## Preserve class and attributes etc.
        })
      }

      # NOTE: It is not sure that the processor supports URLs
      pathnameR <- processor(pathname, ...);

      ## Check if further postprocessoring should be disabled
      metadataR <- getMetadata(pathnameR)
      postprocessR <- getMetadata(pathnameR, "postprocess", local=TRUE)
      if (isFALSE(postprocessR)) metadata$postprocess <- FALSE

      # Always return the relative path
      pathnameR <- getAbsolutePath(pathnameR);
      res <- RspFileProduct(pathnameR, attrs=list(metadata=metadata), mustExist=FALSE);
      res <- setMetadata(res, name="source", value=pathname);

      res;
    } # fcn()
    verbose && cat(verbose, "Processor found: ", type);
  }

  verbose && exit(verbose);

  fcn;
}, protected=TRUE) # findProcessor()



############################################################################
# HISTORY:
# 2015-05-11
# o Added postprocessor for PDF compression.
# 2015-02-04
# o Now the processor function returned by findProcessor() passes
#   meta data as a named list to the underlying compiler/processor.
# 2014-10-18
# o Now the LaTeX processor returned by findProcess() for RspFileProduct
#   (with content type application/x-tex or application/x-latex) will add
#   the directory of the RSP source file to the TEXINPUTS.
# o Now the processor returned by findProcess() for RspFileProduct always
#   set the 'source' metadata.
# 2014-04-18
# o Added argument 'default' to getType() for RspFileProduct.
# o BUG FIX: RspFileProduct would corrupt URLs.
# 2014-03-24
# o WORKAROUND: Due to limitations on browseURL(), view() for
#   RspFileProduct failed to open files with commas in their filenames
#   (e.g. path/foo,bar.html) on some file system (e.g. Windows).  As a
#   workaround, view() now sets the working directory temporarily to that
#   of the file to be displayed and then calls browseURL().
# 2014-02-07
# o Now print() for RspFileProduct reports the file size also in units
#   of kB, MB, GB, etc.
# o Added getFileSize() for RspFileProduct.  Taken from R.filesets.
# 2013-12-14
# o Now getType() for RspFileProduct works also for URLs.
# 2013-03-29
# o Added view().
# 2013-03-25
# o Added Markdown processor.
# 2013-02-18
# o Added argument 'fake' to the returned processor.
# 2013-02-13
# o Added RspProduct and RspFileProduct with corresponding
#   process() methods.
# o Created.
############################################################################
