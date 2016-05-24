###########################################################################/**
# @set class=AbstractCBS
# @RdocMethod report
#
# @title "Generates a report of the segmentation results"
#
# \description{
#  @get "title".
#  Currently reports can be generated for segmentation results of class
#  @see "CBS" and @see "PairedPSCBS".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{An @see "AbstractCBS" object.}
#   \item{sampleName}{A @character string specifying the name of the
#      sample segmented.}
#   \item{studyName}{A @character string specifying the name of study/project.}
#   \item{...}{Optional arguments passed to the RSP template.}
#   \item{rspTags}{Optional @character @vector of tags for further specifying
#      which RSP report to generate.}
#   \item{rootPath}{The root directory where to write the report.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname of the generated PDF.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("report", "AbstractCBS", function(fit, sampleName=getSampleName(fit), studyName, ..., rspTags=NULL, rootPath="reports/", .filename="*", skip=TRUE, envir=new.env(), verbose=FALSE) {
  use("R.rsp (>= 0.20.0)");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  if (is.na(sampleName)) {
    throw("Cannot generate report. Argument 'sampleName' is non-valid or missing.");
  }

  # Argument 'studyName':
  if (missing(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is missing.");
  }
  studyName <- Arguments$getCharacter(studyName);
  if (is.na(studyName)) {
    throw("Cannot generate report. Argument 'studyName' is non-valid.");
  }

  # Argument 'rspTags':
  if (!is.null(rspTags)) {
    rspTags <- Arguments$getCharacters(rspTags);
    rspTags <- unlist(strsplit(rspTags, split=",", fixed=TRUE));
    rspTags <- rspTags[nchar(rspTags) > 0L];
  }

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument '.filename':
  if (!is.null(.filename)) {
    .filename <- Arguments$getCharacter(.filename, useNames=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating CBS report");

  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes(fit));
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments(fit));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Report template arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default arguments
  rspArgs <- list(
    fit = fit,
    sampleName = sampleName,
    studyName = studyName,
    dataSet = NULL,
    Clim = c(0,2*ploidy(fit)),
    Blim = c(0,1),
    figForce = FALSE
  );

  # Override with user arguments
  userArgs <- list(...);
  for (key in names(userArgs)) {
    rspArgs[[key]] <- userArgs[[key]];
  }

  if (is.null(rspArgs$reportPath)) {
    rspArgs$reportPath <- file.path(rootPath, rspArgs$studyName);
  }
  rspArgs$reportPath <- Arguments$getWritablePath(rspArgs$reportPath);
  verbose && cat(verbose, "Report root path: ", rspArgs$reportPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Linking to report files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Linking to report files");

  # Directory where all report templates lives
  srcPath <- "templates";

  # If missing, default to one that comes with PSCBS/templates/
  if (!isDirectory(srcPath)) {
    srcPath <- system.file("templates", package="PSCBS");
  }
  srcPath <- file.path(srcPath, "rsp");
  srcPath <- Arguments$getReadablePath(srcPath);
  verbose && cat(verbose, "Source path: ", srcPath);

  filename <- .filename;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create file links to the main RSP report template
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Main RSP template");

  # Construct the filename of the main RSP file to compile, iff missing.
  if (filename == "*") {
    className <- class(fit)[1];
    fullname <- paste(c(className, rspTags), collapse=",");
    filename <- sprintf("%s,report.tex.rsp", fullname);
  }

  rspPathname <- file.path(srcPath, filename);
  verbose && cat(verbose, "RSP report template: ", rspPathname);
  rspPathname <- Arguments$getReadablePathname(rspPathname);

  destFilename <- sprintf("%s,%s", sampleName, filename);
  destPathname <- filePath(rspArgs$reportPath, destFilename);
  target <- rspPathname;
  link <- destPathname;
  if (!isFile(link)) {
    verbose && cat(verbose, "Adding link: ", link, " -> ", target);
    createLink(link=link, target=target);
  }
  # Sanity check
  stopifnot(isFile(link));

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Skip?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (skip) {
    ## Try to guess report filename
    filename <- basename(destPathname)
    filename <- gsub("[.]rsp$", "", filename)
    ext <- gsub(".*[.]", "", filename)
    fullname <- gsub(sprintf("[.]%s$", ext), "", filename)
    ext <- switch(ext, tex="pdf", md="html", ext)
    filename <- sprintf("%s.%s", fullname, ext)
    pathname <- file.path(rspArgs$reportPath, filename)
    pathname <- getAbsolutePath(pathname)
    verbose && cat(verbose, "Expected output pathname: ", pathname);
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already exists: Skipping.")
      report <- R.rsp::RspFileProduct(pathname)
      return(report)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create file links to all LaTeX include files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "All LaTeX files");

  files <- list.files(path=srcPath, pattern="[.](bib|bst|cls|sty|tex)$", full.names=TRUE, recursive=FALSE);
  files <- files[file_test("-f", files)];
  if (length(files) > 0L) {
    verbose && cat(verbose, "Number of such files found: ", length(files));
    verbose && print(verbose, files);

    for (kk in seq_along(files)) {
      target <- files[kk];
      link <- filePath(rspArgs$reportPath, basename(files[kk]));
      if (!isFile(link)) {
        verbose && cat(verbose, "Adding link: ", link, " -> ", target);
        createLink(link=link, target=target);
      }
      # Sanity check
      stopifnot(isFile(link));
    }
  } else {
    verbose && cat(verbose, "No such files found.");
  }

  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create file links to all 'incl.*' subdirectories
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "All 'incl.*' subdirectories");

  dirs <- list.files(srcPath, pattern="^incl", full.names=TRUE, recursive=FALSE);
  dirs <- dirs[file_test("-d", dirs)];

  if (length(dirs) > 0L) {
    verbose && cat(verbose, "Number of such directories found: ", length(dirs));
    verbose && print(verbose, dirs);

    for (kk in seq_along(dirs)) {
      target <- dirs[kk];
      link <- filePath(rspArgs$reportPath, basename(dirs[kk]));
      if (!isDirectory(link)) {
        verbose && cat(verbose, "Adding link: ", link, " -> ", target);
        createLink(link=link, target=target);
      }
      # Sanity check
      stopifnot(isDirectory(link));
    }
  } else {
    verbose && cat(verbose, "No such directories found.");
  }

  verbose && exit(verbose);


  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build reports
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Processing RSP template");
  rspArgs$figPath <- "figures/";
  args <- c(list(rspArgs=rspArgs), rspArgs);
  report <- R.rsp::rfile(destPathname, workdir=rspArgs$reportPath, args=args, envir=envir, verbose=verbose);
  verbose && exit(verbose);

  verbose && cat(verbose, "Final report: ", getRelativePath(report));

  verbose && exit(verbose);

  report;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-03-21
# o Now report() uses file links instead of copying template files.
#   It also links to all LaTeX related files and all directories
#   named '^incl.*'.
# 2013-03-09
# o Now report() also included files listed in the optional file
#   '.install_extras' of the source RSP template directory.
#   The same filename is used by 'R CMD build/check' for including
#   extra vignette source files.
# 2012-09-18
# o Added argument 'force' to report() for AbstractCBS.  This will
#   copy the RSP template files again, although they are already in
#   reports/ output directory.
# o Now report(fit, ..., rspTags) for AbstractCBS looks for the RSP
#   template named <className>(,<rspTags>),report.tex.rsp, where
#   className is class(fit)[1] and  argument 'rspTags' is an optional
#   comma-separated character string/vector.
# o Now report() for AbstractCBS looks for the RSP template in templates/,
#   and as a backup in templates,PSCBS/.  If the latter does not exist,
#   it is automatically created as a soft link to templates/ of the
#   PSCBS package.  This allows anyone to create their own customized
#   copy (in templates/) of the default PSCBS RSP report.
# 2012-05-30
# o Now report() gives more a informative error message if arguments
#   'sampleName' or 'studyName' are non-valid or missing.
# 2012-04-20
# o Added argument '.filenames'.
# o Created from former PairedPSCBS.REPORT.R, which history as below.
# 2012-02-27
# o Added Rdoc help.
# o Added report() for PairedPSCBS.
# o Created.
############################################################################
