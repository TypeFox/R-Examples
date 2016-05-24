###########################################################################/**
# @RdocDefault rfile
# @alias rfile.RspString
# @alias rfile.RspDocument
# @alias rfile.RspRSourceCode
# @alias rfile.function
# @alias rfile.expression
#
# @title "Evaluates and postprocesses an RSP document and outputs the final RSP document file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{file, path}{Specifies the RSP file to processed, which can
#      be a file, a URL or a @connection.
#      If a file, the \code{path} is prepended to the file, iff given.}
#   \item{output}{A @character string or a @connection specifying where
#      output should be directed.
#      The default is a file with a filename where the file extension
#      (typically \code{".rsp"}) has been dropped from \code{file}
#      in the directory given by the \code{workdir} argument.}
#   \item{workdir}{The working directory to use after parsing and
#      preprocessing, but while \emph{evaluating} and \emph{postprocessing}
#      the RSP document.
#      If argument \code{output} specifies an absolute pathname,
#      then the directory of \code{output} is used, otherwise the
#      current directory is used.}
#   \item{type}{The default content type of the RSP document.  By default, it
#      is inferred from the \code{output} filename extension, iff possible.}
#   \item{envir}{The @environment in which the RSP document is
#      preprocessed and evaluated.}
#   \item{args}{A named @list of arguments assigned to the environment
#     in which the RSP string is parsed and evaluated.
#     See @see "R.utils::cmdArgs".}
#   \item{postprocess}{If @TRUE, and a postprocessing method exists for
#      the generated RSP product, it is postprocessed as well.}
#   \item{...}{Additional arguments passed to the RSP engine.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an @see "RspProduct".
#   If argument \code{output} specifies a file, then this is
#   an @see "RspFileProduct".
# }
#
# \section{Processing RSP files from the command line}{
#   Using @see "Rscript" and \code{rfile()}, it is possible to process
#   an RSP file from the command line.  For example,
#
#   \code{Rscript -e "R.rsp::rfile(file='RSP_refcard.tex.rsp', path=system.file('doc', package='R.rsp'))"}
#
#   parses and evaluates \file{RSP_refcard.tex.rsp} and output \file{RSP_refcard.pdf} in the current directory.
# }
#
# \examples{
# @include "../incl/rfile.Rex"
#
# \dontrun{
# # Compile and display the main vignette (requires LaTeX)
# if (isCapableOf(R.rsp, "latex")) {
#   path <- system.file("doc", package="R.rsp")
#   pdf <- rfile("Dynamic_document_creation_using_RSP.tex.rsp", path=path)
#   cat("Created document: ", pdf, "\n", sep="")
#   if (interactive()) browseURL(pdf)
# }
# }
# }
#
# @author
#
# \seealso{
#  @see "rstring" and @see "rcat".
# }
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("rfile", "default", function(file, path=NULL, output=NULL, workdir=NULL, type=NA, envir=parent.frame(), args="*", postprocess=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file' & 'path':
  if (inherits(file, "connection")) {
  } else if (inherits(file, "RspFileProduct")) {
  } else if (is.character(file)) {
    if (!is.null(path)) {
      file <- file.path(path, file);
    }
    if (!isUrl(file)) {
      withoutGString({
        file <- Arguments$getReadablePathname(file, absolute=TRUE);
      })
    }
  }

  # Argument 'workdir':
  if (is.null(workdir)) {
    if (isAbsolutePath(output)) {
      workdir <- getParent(output);
    } else {
      workdir <- ".";
    }
  }
  workdir <- Arguments$getWritablePath(workdir);
  if (is.null(workdir)) workdir <- ".";

  # Argument 'output':
  if (is.null(output)) {
    if (inherits(file, "connection")) {
      throw("When argument 'file' is a connection, then 'output' must be specified.");
    }

    # Is the input a filename or an URI?
    if (isUrl(file)) {
      # If URI, drop any URI arguments
      url <- splitUrl(file);
      filename <- basename(url$path);
      filename <- Arguments$getReadablePathname(filename, adjust="url", mustExist=FALSE);
    } else {
      filename <- basename(file);
    }

    pattern <- "((.*)[.]([^.]+)|([^.]+))[.]([^.]+)$";
    outputF <- gsub(pattern, "\\1", filename, ignore.case=TRUE);
    withoutGString({
      output <- Arguments$getWritablePathname(outputF, path=workdir);
    })
    output <- getAbsolutePath(output);
    # Don't overwrite the input file
    if (output == file) {
      throw("Cannot process RSP file. The inferred argument 'output' is the same as argument 'file' & 'path': ", output, " == ", file);
    }
  } else if (inherits(output, "connection")) {
  } else if (identical(output, "")) {
    output <- stdout();
  } else if (is.character(output)) {
    withoutGString({
      if (isAbsolutePath(output)) {
        output <- Arguments$getWritablePathname(output);
      } else {
        output <- Arguments$getWritablePathname(output, path=workdir);
        output <- getAbsolutePath(output);
      }
    })
    if (is.character(file) && (output == file)) {
      throw("Cannot process RSP file. Argument 'output' specifies the same file as argument 'file' & 'path': ", output, " == ", file);
    }
  } else {
    throw("Argument 'output' of unknown type: ", class(output)[1L]);
  }

  # Argument 'type':
  if (is.null(type)) {
    if (is.character(output)) {
      type <- extensionToIMT(output);
      attr(type, "fixed") <- TRUE;
    } else {
      type <- NA;
    }
  }
  if (is.na(type)) {
    if (is.character(output)) {
      type <- extensionToIMT(output);
    }
  }
  fixed <- attr(type, "fixed");
  type <- Arguments$getCharacter(type);
  attr(type, "fixed") <- fixed;

  # Argument 'envir':
  stopifnot(is.environment(envir));

  # Argument 'args':
  args <- cmdArgs(args=args);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Processing RSP file");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Information on input and output
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (verbose) {
    # Information on input
    if (inherits(file, "RspFileProduct")) {
      cat(verbose, "Input file:");
      print(verbose, file);
    } else if (is.character(file)) {
      if (isUrl(file)) {
        cat(verbose, "Input URL: ", file);
      } else {
        cat(verbose, "Input pathname: ", file);
      }
    } else if (inherits(file, "connection")) {
      ci <- summary(file);
      printf(verbose, "Input '%s' connection: %s\n",
          class(ci)[1L], ci$description);
    }

    # Information on output
    if (is.character(output)) {
      cat(verbose, "Output pathname: ", output);
    } else if (inherits(output, "connection")) {
      ci <- summary(output);
      printf(verbose, "Output '%s' connection: %s\n",
          class(ci)[1L], ci$description);
    }

    # Information on content *output* type
    printf(verbose, "Default content type: %s\n", type);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assign RSP arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Assigning RSP arguments");
  verbose && cat(verbose, "Environment: ", getName(envir));
  if (length(args) > 0L) {
    verbose && cat(verbose, "Arguments assigned: ", hpaste(names(args)));
    # Assign arguments to the parse/evaluation environment
    attachLocally(args, envir=envir);
  } else {
    verbose && cat(verbose, "Arguments assigned: <none>");
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce to an RspFileProduct
  # FIXME: Would be nice to be able to handle unknown file extensions,
  # e.g. when an RSP file is downloaded from online wihtout a filename ext.
  if (!inherits(file, "RspFileProduct")) {
    file <- RspFileProduct(file, mustExist=FALSE);
    processor <- findProcessor(file);
    if (is.null(processor)) {
      # Assume an RSP document if type cannot be inferred by filename etc.
      file <- RspFileProduct(file, type="application/x-rsp", mustExist=FALSE);
    }
  }

  # Process...
  if (getType(file, default="text/plain") == "application/x-rsp") {
    # (a) An RSP document, or...
    verbose && enter(verbose, "Reading RSP document");
    str <- .readText(file);
    verbose && printf(verbose, "Number of characters: %d\n", nchar(str));
    verbose && str(verbose, str);
    verbose && exit(verbose);

    verbose && enter(verbose, "Parsing RSP document");
    rstr <- RspString(str, type=type, source=file);
    rstr <- setMetadata(rstr, name="source", value=file);
    doc <- parse(rstr, envir=envir, ...);
    verbose && print(verbose, doc);
    rstr <- str <- NULL; # Not needed anymore
    verbose && exit(verbose);

    res <- rfile(doc, output=output, workdir=workdir, envir=envir, args=NULL, postprocess=postprocess, ..., verbose=verbose);
  } else {
    # (b) ...other type of document.
    res <- process(file, workdir=workdir, envir=envir, args=NULL, recursive=postprocess, ..., verbose=verbose);
    res <- setMetadata(res, name="source", value=file);
  }

  verbose && exit(verbose);

  res;
}) # rfile()



setMethodS3("rfile", "RspString", function(rstr, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Processing RSP string");

  verbose && enter(verbose, "Parsing RSP document");
  doc <- parse(rstr, ...);
  verbose && print(verbose, doc);
  rstr <- str <- NULL; # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Translating RSP document (to R)");
  rcode <- toR(doc, ...);
  verbose && printf(verbose, "Number of R source code lines: %d\n", length(rcode));
  doc <- NULL; # Not needed anymore
  verbose && exit(verbose);

  res <- rfile(rcode, ..., verbose=verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # rfile()



setMethodS3("rfile", "RspDocument", function(doc, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Processing RSP document");

  verbose && enter(verbose, "Translating RSP document (to R)");
  rcode <- toR(doc, ...);
  verbose && printf(verbose, "Number of R source code lines: %d\n", length(rcode));
  doc <- NULL; # Not needed anymore
  verbose && exit(verbose);

  res <- rfile(rcode, ..., verbose=verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # rfile()


setMethodS3("rfile", "RspRSourceCode", function(rcode, output, workdir=NULL, envir=parent.frame(), args="*", postprocess=TRUE, ..., verbose=FALSE) {
  # In-string variable substitute
  vsub <- function(pathname, ...) {
    gstr <- GString(pathname);
    str <- gstring(gstr, where=c("envir", "Sys.getenv", "getOption"),
                         envir=envir, inherits=FALSE)[1L];
    str <- wstring(str, envir=envir);
    str;
  } # vsub()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'workdir':
  if (is.null(workdir)) {
    if (isAbsolutePath(output)) {
      workdir <- getParent(output);
    } else {
      workdir <- ".";
    }
  }
  workdir <- Arguments$getWritablePath(workdir);
  if (is.null(workdir)) workdir <- ".";

  # Argument 'output':
  if (inherits(output, "connection")) {
  } else if (identical(output, "")) {
    output <- stdout();
  } else if (is.character(output)) {
    withoutGString({
      if (isAbsolutePath(output)) {
        output <- Arguments$getWritablePathname(output);
      } else {
        output <- Arguments$getWritablePathname(output, path=workdir);
        output <- getAbsolutePath(output);
      }
    })
  } else {
    throw("Argument 'output' of unknown type: ", class(output)[1L]);
  }

  # Argument 'envir':
  stopifnot(is.environment(envir));

  # Argument 'args':
  args <- cmdArgs(args=args);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Processing RSP R source code");

  if (verbose) {
    if (is.character(output)) {
      cat(verbose, "Output pathname: ", output);
    } else if (inherits(output, "connection")) {
      ci <- summary(output);
      printf(verbose, "Output '%s' connection: %s\n",
          class(ci)[1L], ci$description);
    }
  }


  verbose && enter(verbose, "Assigning RSP arguments");
  verbose && cat(verbose, "Environment: ", getName(envir));
  if (length(args) > 0L) {
    verbose && cat(verbose, "Arguments assigned: ", hpaste(names(args)));
    # Assign arguments to the parse/evaluation environment
    attachLocally(args, envir=envir);
  } else {
    verbose && cat(verbose, "Arguments assigned: <none>");
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Evaluating RSP R source code");

  # Change working directory
  opwd <- NULL;
  if ((workdir != ".") && (workdir != getwd())) {
    opwd <- getwd();
    on.exit({ if (!is.null(opwd)) setwd(opwd) }, add=TRUE);
    verbose && cat(verbose, "Temporary working directory: ", getAbsolutePath(workdir));
    setwd(workdir);
  }

  res <- rcat(rcode, output=output, envir=envir, args=NULL, ..., verbose=less(verbose, 10));

  withoutGString({
    if (isFile(output)) {
      res <- RspFileProduct(output, attrs=getAttributes(res));

      # Rename output file via GString substitution of the filename?
      resG <- vsub(res);
      if (resG != res) {
        if (renameFile(res, resG, overwrite=TRUE)) {
          # FIXME: res <- newInstance(res, resG);
          res <- RspFileProduct(resG, attrs=getAttributes(res));
        } else {
          warning(sprintf("Failed to rename output file containing variable substitutions in its name (keeping the current one): ", sQuote(res), " -> ", sQuote(resG)));
        }
      }
      resG <- NULL; # Not needed anymore
    } else {
      res <- RspProduct(output, attrs=getAttributes(res));
    }
  }) # withoutGString()

  verbose && print(verbose, res);
  rcode <- output <- NULL; # Not needed anymore

  # Reset the working directory?
  if (!is.null(opwd)) {
    setwd(opwd);
    opwd <- NULL;
  }

  verbose && exit(verbose);

  if (postprocess && hasProcessor(res)) {
    res <- process(res, workdir=workdir, ..., verbose=verbose);
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE) # rfile()


setMethodS3("rfile", "function", function(object, ..., envir=parent.frame(), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'object':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rfile() for ", class(object)[1L]);

  ## Temporarily assign the function to the evaluation environment
  ## and set its own environment also to the evaluation environment
  fcn <- object;
  environment(fcn) <- envir;
  fcnName <- tempvar(".rfcn", value=fcn, envir=envir);
  on.exit({
    rm(list=fcnName, envir=envir, inherits=FALSE);
  }, add=TRUE);
  code <- sprintf("%s()", fcnName);
  rcode <- RspRSourceCode(code);
  res <- rfile(rcode, ..., envir=envir, verbose=verbose);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # rfile()


setMethodS3("rfile", "expression", function(object, ..., envir=parent.frame(), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'object':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "rfile() for ", class(object)[1L]);
  # Deparsing 'object[[1L]]' instead of 'object' in order to drop
  # the 'expression({ ... })' wrapper.
  code <- deparse(object[[1L]]);
  rcode <- RspRSourceCode(code);
  res <- rfile(rcode, ..., envir=envir, verbose=verbose);
  verbose && exit(verbose);

  res;
}, protected=TRUE) # rfile()


############################################################################
# HISTORY:
# 2015-02-14
# o CRAN: Used to have \donttest{} in example(rfile), which was there to
#   avoid the test running longer than 5 seconds.  This was disapproved
#   in resubmission.  Now using \dontrun{} instead.
# 2014-05-30
# o Now metadata 'source' is set by rfile(), iff possible.  It gives the
#   absolute path to the input file, or the URL, of the source RSP file.
#   It can be accessed via <%@meta name="source"%>.
# 2014-01-02
# o Added rstring(), rcat() and rfile() for expression:s too.
# 2013-12-14
# o Now rfile() accepts also non-RSP documents, e.g. rfile("report.md"),
#   rfile("report.Rnw"), and rfile("report.tex").
# 2013-12-13
# o Now rfile() does string variable substitution of the output pathname,
#   if possible.
# o Now rfile() no longer interprets input pathnames as GString:s.
# 2013-07-16
# o Added rstring(), rcat() and rfile() for function:s.
# 2013-07-14
# o Added rfile() for RspSourceCode, which now is utilized by the default
#   rfile().
# 2013-05-22
# o ROBUSTNESS: Now rfile() handles files with only one filename extension.
# 2013-02-23
# o Now rfile() can also infer default filenames from URLs with parameters.
# 2013-02-18
# o Added argument 'fake' to rfile().
# 2013-02-13
# o Added argument 'postprocess' to rfile().
# o Now rfile() sets the default content type based on the filename
#   extension, iff possible.
# o Added argument 'workdir' to rfile().
# o Added support for 'pathname' also being a connection or a URL.
# o Added some protection against overwriting the input file.
# o Renamed rspPlain() to rfile(), cf. rstring() and rcat().
# o rspPlain() is now utilizing the new RSP engine, e.g. rcat().
# o CLEANUP: Removed all dependencies on RspResponse, FileRspResponse etc.
# 2013-02-08
# o Extracted rspPlain() from rsp().
# o Created.
############################################################################
