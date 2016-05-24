###########################################################################/**
# @set class=GenericDataFileSet
# @RdocMethod dsApply
# @aliasmethod dsApplyInPairs
# @alias dsApplyInPairs
#
# @title "Applies a function to each file in the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{ds, ds1, ds2}{@see "GenericDataFileSet":s.}
#  \item{IDXS}{A (named) @list, where each element contains a @vector data
#    set indices, or an @integer @vector of individual elements.
#    If @NULL, then ... with names as of the data set.}
#  \item{DROP}{If @FALSE, the first argument passed to \code{FUN} is always a @list of files.
#    If @TRUE, an single-index element is passed to \code{FUN} as a file instead of
#    as a @list containing a single file.}
#  \item{AS}{(optional) A @function coercing the first set/group object passed.}
#  \item{FUN}{A @function.}
#  \item{...}{Arguments passed to \code{FUN}.}
#  \item{args}{(optional) A @list of additional arguments
#    passed to \code{FUN}.}
#  \item{skip}{If @TRUE, already processed files are skipped.}
#  \item{verbose}{See @see "R.utils::Verbose".}
#  \item{.parallel}{A @character string specifying what mechanism to use
#    for performing parallel processing, if at all.}
#  \item{.control}{(internal) A named @list structure controlling
#        the processing.}
# }
#
# \value{
#   Returns a named @list where the names are those of argument \code{IDXS}.
# }
#
# \examples{\dontrun{
#  @include "../incl/GenericDataFileSet.dsApply.Rex"
# }}
#
# \seealso{
#   The \pkg{future}, \pkg{BiocParallel} and \pkg{BatchJobs} packages
#   are utilized for parallel/distributed processing, depending on settings.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("dsApply", "GenericDataFileSet", function(ds, IDXS=NULL, DROP=is.null(IDXS), AS=as.list, FUN, ..., args=list(), skip=FALSE, verbose=FALSE, .parallel=c("none", "future", "BatchJobs", "BiocParallel::BatchJobs"), .control=list(dW=1.0)) {
  # To please R CMD check, because
  # (i) BatchJobs is just "suggested"
  getJobNr <- batchMap <- showStatus <- findNotSubmitted <-
      findNotRunning <- submitJobs <- findNotTerminated <-
      loadResults <- loadConfig <- NULL;
  # (ii) BiocParallel is just "suggested"
  BatchJobsParam <- register <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  assertNoGlobalVariables <- function(FUN, ...) {
    # TO DO...
    ## globals <- findGlobals(FUN, merge=FALSE);
  } # assertNoGlobalVariables()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'DROP':
  DROP <- Arguments$getLogical(DROP);

  # Argument 'IDXS':
  if (is.null(IDXS)) {
    IDXS <- seq_along(ds);
    names(IDXS) <- getFullNames(ds);
    IDXS <- as.list(IDXS);
  } else if (is.numeric(IDXS)) {
    max <- length(ds);
    IDXS <- Arguments$getIndices(IDXS, max=max);
    if (is.null(names(IDXS))) {
      names(IDXS) <- getFullNames(ds);
      IDXS <- as.list(IDXS);
    }
  } else if (is.list(IDXS)) {
    max <- length(ds);
    for (idxs in IDXS) {
      idxs <- Arguments$getIndices(idxs, max=max);
    }
  } else {
    throw("Invalid argument 'IDXS': ", class(IDXS)[1L]);
  }

  # Argument 'FUN':
  stopifnot(is.function(FUN));
  assertNoGlobalVariables(FUN);


  # Arguments '...':
  vargs <- list(...);
  nvargs <- length(vargs);

  # Argument 'args':
  if (!is.list(args)) {
    throw("Argument 'args' must be a list: ", mode(args));
  }

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument '.parallel':
  if (missing(.parallel)) {
    .parallel <- getOption("R.filesets/parallel", "none")
  }
  parallel <- match.arg(.parallel)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Processing ", class(ds)[1L]);
  verbose && cat(verbose, "Mechanism for parallel processing: ", parallel);
  verbose && print(verbose, ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sets of files to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Number of subsets to be processed: ", length(IDXS));
  verbose && str(verbose, head(IDXS));

  # Drop only if all index groups have exactly one file
  if (DROP) {
    ns <- sapply(IDXS, FUN=length);
    DROP <- all(ns == 1L);
  }

  # Still drop?
  if (DROP) {
    sets <- lapply(IDXS, FUN=function(idx) ds[[idx]]);
  } else {
    sets <- vector("list", length=length(IDXS));
    for (gg in seq_along(IDXS)) {
      idxs <- IDXS[[gg]];
      set <- ds[idxs];
      # FIXME/BACKWARD COMPATIBLE? /HB 2014-03-30
      if (is.function(AS)) {
        set <- AS(set);
        if (identical(AS, as.list)) {
          name <- names(IDXS)[gg];
          if (!is.null(name)) names(set)[1L] <- name;
        }
      }
      sets[[gg]] <- set;
    } # for (gg ...)
  }
  names(sets) <- names(IDXS);
  if (is.null(names(sets))) {
    names(sets) <- sprintf("<Group #%d>", seq_along(sets));
  }
  verbose && str(verbose, head(sets));
  IDXS <- NULL; # Not needed anymore

  # FIXME/AD HOC/BACKWARD COMPATIBLE /HB 2014-03-30
  # Set attribute 'groupName' for each element.  This is used to pass
  # the group name to FUN().
  for (gg in seq_along(sets)) {
    set <- sets[[gg]];
    name <- names(sets)[gg];
    attr(set, "name") <- name;
    sets[[gg]] <- set;
  }

  # The additional set of arguments passed in each function call
  vargs <- c(vargs, args);
  allArgs <- c(vargs, list(skip=skip, verbose=verbose));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 1: Sequentially using regular R
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- NULL;
  if (parallel == "none") {
    # Allocate result list
    res <- vector("list", length=length(sets));
    names(res) <- names(sets);

    for (gg in seq_along(sets)) {
      name <- names(sets)[gg];
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d", gg, name, length(sets)));
      set <- sets[[gg]];
      verbose && print(verbose, set);
      argsGG <- c(list(set), allArgs);
      verbose && cat(verbose, "Call arguments:");
      argsT <- argsGG; argsT$verbose <- as.character(argsT$verbose);
      verbose && str(verbose, argsT);
      argsT <- NULL; # Not needed anymore
      resGG <- do.call(FUN, args=argsGG);
      verbose && str(verbose, resGG);

      # Record
      res[[gg]] <- resGG;

      # Not needed anymore
      idxs <- set <- argsGG <- resGG <- NULL;

      verbose && exit(verbose);
    } # for (gg ...)
  } # if (parallel == "none")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 2: Evaluation using futures
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel == "future") {
    verbose && enter(verbose, "Processing using futures");

    # Allocate result list
    futures <- listenv()

    for (gg in seq_along(sets)) {
      name <- names(sets)[gg]
      verbose && enter(verbose, sprintf("Group #%d ('%s') of %d", gg, name, length(sets)))
      set <- sets[[gg]]
      verbose && print(verbose, set)
      argsGG <- c(list(set), allArgs)
      verbose && cat(verbose, "Call arguments:")
      argsT <- argsGG; argsT$verbose <- as.character(argsT$verbose)
      verbose && str(verbose, argsT)
      argsT <- NULL; # Not needed anymore

      futureGG <- future(do.call(FUN, args=argsGG))
      verbose && str(verbose, futureGG)

      # Record
      futures[[gg]] <- futureGG

      # Not needed anymore
      idxs <- set <- argsGG <- futureGG <- NULL

      verbose && exit(verbose)
    } # for (gg ...)

    ## No longer needed
    rm(list=c("FUN", "argsGG"))

    names(futures) <- names(sets)

    ## Resolve the value of all futures
    res <- lapply(futures, FUN=value)

    ## Not needed anymore
    rm(list="futures")
    verbose && exit(verbose)
  } # if (parallel == "future")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 3: BatchJobs processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel == "BatchJobs") {
    verbose && enter(verbose, "Processing using BatchJobs");

    # Attach "suggested" BatchJobs package
    .useBatchJobs();

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setup
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Poll roughly every dW seconds
    dW <- .control$dW;
    if (is.null(dW)) dW <- 1.00;
    dW <- Arguments$getNumeric(dW, range=c(0,Inf));

    # BatchJob registry to be used
    reg <- .getBatchJobRegistry(ds, args=vargs);
    verbose && print(verbose, reg);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (i) Add jobs, iff missing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nbrOfJobs <- getJobNr(reg);
    if (nbrOfJobs == 0L) {
      verbose && enter(verbose, "Adding jobs to registry");
      # Tweak FUN()
##      FUNx <- function(...) {
##        # Allow comma-separated registry directory paths
##        oopts <- options("BatchJobs.check.posix");
##        on.exit({ options(oopts) }, add=TRUE);
##        options("BatchJobs.check.posix"=FALSE);
##        FUN(...);
##      } # FUNx()
      ids <- batchMap(reg, fun=FUN, sets, more.args=allArgs);
      verbose && cat(verbose, "Job IDs added:");
      verbose && str(verbose, ids);
      verbose && print(verbose, reg);
      verbose && exit(verbose);
    }

    verbose && print(verbose, showStatus(reg));
#    throw("Jobs have already been added: ", reg$id);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (ii) Launch jobs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Launching jobs");
    lastTodo <- NULL;
    todo <- findNotSubmitted(reg);
    if (length(todo) > 0L) {
      # (a) Wait and see if jobs are being submitted by other process
      while(!identical(todo, lastTodo)) {
         lastTodo <- todo;
         Sys.sleep(1.0);
         todo <- findNotRunning(reg);
      }

      verbose && cat(verbose, "Job IDs to be submitted:");
      verbose && print(verbose, todo);
      submitted <- submitJobs(reg, ids=todo);
      verbose && cat(verbose, "Job IDs actually submitted:");
      verbose && print(verbose, submitted);
      verbose && cat(verbose, "Job IDs not submitted:");
      verbose && print(verbose, setdiff(todo, submitted));
    } else {
      verbose && cat(verbose, "No new jobs to be submitted.");
    }

    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (iii) Wait for jobs to finish
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Waiting for jobs to finish");
    t0 <- Sys.time();
    tCount <- 0L;
    status <- NULL;
    while (length(findNotTerminated(reg)) > 0L) {
      lastStatus <- status;
      out <- capture.output(status <- showStatus(reg));
      if (identical(status, lastStatus)) {
        verbose && writeRaw(verbose, ".");
        # Time stamp?
        dt <- difftime(Sys.time(), t0, units="secs");
        dMins <- as.integer(dt) %/% 10;
        if (dMins > tCount) {
          tCount <- dMins;
          if (dt > 1.5*60) {
            units(dt) <- "mins";
          } else if (dt > 1.5*3600) {
            units(dt) <- "hours";
          }
          verbose && writeRaw(verbose, sprintf("[%s]\n", format(dt)));
        }
      } else {
        verbose && writeRaw(verbose, "\n");
        verbose && print(verbose, status);
      }
      Sys.sleep(dW);
    } # while(...)
    verbose && exit(verbose);

    verbose && print(verbose, showStatus(reg));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (iv) Retrieve results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Retrieve results");
    res <- loadResults(reg, simplify=FALSE, missing.ok=TRUE);
    names(res) <- names(sets);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (parallel == "BatchJobs")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 4: BiocParallel w/ BatchJobs processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel == "BiocParallel::BatchJobs") {
    verbose && enter(verbose, "Processing using BiocParallel");

    # WORKAROUND: Make sure 'methods' package is *attached*, not
    # just loaded. /HB 2013-11-09
    .require <- require   # To please R CMD check
    .require("methods") || throw("Package not attached: methods")

    # Attach "suggested" BiocParallel package
    .useBiocParallel();

    # Attach "suggested" BatchJobs package
    .useBatchJobs();

    conffile <- c(".BatchJobs.R", "~/.BatchJobs.R")
    conffile <- normalizePath(conffile, mustWork=FALSE);
    conffile <- conffile[file_test("-f", conffile)];
    if (length(conffile) > 0L) {
      conffile <- conffile[1L];
      if (isFile(conffile)) loadConfig(conffile);
    }

    # WORKAROUND: Allow for commas in BatchJobs-related pathnames
    oopts <- options("BatchJobs.check.posix");
    on.exit({ options(oopts) }, add=TRUE);
    options("BatchJobs.check.posix"=FALSE);

    # It's here one can specify PBS options such as number of
    # nodes, number of cores, walltime etc.
    bpParam <- BatchJobsParam(resources=NULL);
    register(bpParam);
    verbose && cat(verbose, "Using parameters:");
    verbose && print(verbose, bpParam);

    verbose && enter(verbose, "Calling bplapply()");
    args <- c(list(sets, FUN=FUN), allArgs, BPPARAM=bpParam);
    verbose && cat(verbose, "Arguments passed to bplapply():");
    argsT <- args; argsT$verbose <- as.character(argsT$verbose);
    verbose && str(verbose, argsT);
    argsT <- NULL; # Not needed anymore
    res <- do.call(BiocParallel::bplapply, args);
    names(res) <- names(sets);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (parallel == "BiocParallel")

  verbose && exit(verbose);

  res;
}, protected=TRUE) # dsApply()



setMethodS3("dsApplyInPairs", "GenericDataFileSet", function(ds1, ds2, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ds1' & 'ds2':
  ds2 <- Arguments$getInstanceOf(ds2, class(ds1)[1L])
  stopifnot(length(ds2) == length(ds1))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # dsApply() in pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsP <- c(ds1, ds2)
  idxs <- matrix(seq_along(dsP), nrow=2L, byrow=TRUE)
  IDXS <- as.list(as.data.frame(idxs))
  names <- getFullNames(dsP)
  names(IDXS) <- sapply(IDXS, FUN=function(idxs) {
    sprintf("Pair (%s,%s)", names[idxs[1]], names[idxs[2]])
  })
  idxs <- names <- ds1 <- ds2 <- NULL  # Not needed anymore

  dsApply(dsP, IDXS=IDXS, ...)
}, protected=TRUE) # dsApplyInPairs()



setMethodS3(".getBatchJobRegistryId", "default", function(class, label=NULL, version=NULL, ..., verbose=FALSE) {
  # Argument 'class':
  class <- Arguments$getCharacters(class);

  # Argument 'label':
  label <- Arguments$getCharacter(label);

  # Argument 'version':
  version <- Arguments$getCharacter(version);


  # Construct key from all object/arguments
  key <- list(
    class=class,
    label=label,
    version=version,
    ...
  );
  keyId <- R.cache::getChecksum(key, algo="crc32");
  keyTime <- format(Sys.time(), "%Y%m%d%H%M%S");
  pid <- Sys.getpid();
  id <- paste(c(key$class[1L], key$label, keyId, keyTime, pid), collapse="_");
  id <- Arguments$getCharacter(id, length=c(1L,1L));

  id;
}, private=TRUE)


setMethodS3(".getBatchJobRegistryId", "GenericDataFileSet", function(object, ...) {
  keys <- lapply(object, FUN=function(file) {
    list(filename=getFilename(file), fileSize=getFileSize(file));
  });
  .getBatchJobRegistryId(class=class(object), fileKeys=keys, ...);
}, private=TRUE)



setMethodS3(".getBatchJobRegistry", "default", function(..., skip=TRUE) {
  .useBatchJobs();

  # To please R CMD check (already loaded above)
  requireNamespace("BatchJobs");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'skip':
  skip <- Arguments$getLogical(skip);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constructor BatchJobs registry ID...
  id <- .getBatchJobRegistryId(...);

  # ...and registry directory
  rootPath <- ".batchJobsRegistries";
  path <- file.path(rootPath, id);
  path <- Arguments$getWritablePath(path);

  # Allow comma-separated registry directory paths
  oopts <- options("BatchJobs.check.posix");
  on.exit({ options(oopts) }, add=TRUE);
  options("BatchJobs.check.posix"=FALSE);

  # Load BatchJobs registry or create if missing
  packages <- c("R.filesets");
  reg <- BatchJobs::makeRegistry(id=id, file.dir=path, skip=skip, packages=packages);

  reg;
}, private=TRUE) # .getBatchJobRegistry()


.useBatchJobs <- function(version="*", ...) {
  if (identical(version, "*")) {
    desc <- packageDescription("R.filesets");
    desc <- desc[c("Depends", "Imports", "Suggests")];
    desc <- unlist(desc, use.names=FALSE);
    desc <- strsplit(desc, split=",", fixed=TRUE);
    desc <- unlist(desc, use.names=FALSE);
    desc <- gsub("\n", "", desc, fixed=TRUE);
    pattern <- "[ ]*([^ (]+)[ ]*(()|[(][^]]+[)])";
    pkgs <- gsub(pattern, "\\1", desc);
    vers <- gsub(pattern, "\\2", desc);
    vers <- gsub("[()]", "", vers);
    names(vers) <- pkgs;
    dups <- duplicated(pkgs);
    vers <- vers[!dups];
    version <- vers["BatchJobs"];
  }
  R.utils::use("BatchJobs", version=version, ...);
} # .useBatchJobs()


.useBiocParallel <- function(version="*", ...) {
  if (identical(version, "*")) {
    desc <- packageDescription("R.filesets");
    desc <- desc[c("Depends", "Imports", "Suggests")];
    desc <- unlist(desc, use.names=FALSE);
    desc <- strsplit(desc, split=",", fixed=TRUE);
    desc <- unlist(desc, use.names=FALSE);
    desc <- gsub("\n", "", desc, fixed=TRUE);
    pattern <- "[ ]*([^ (]+)[ ]*(()|[(][^]]+[)])";
    pkgs <- gsub(pattern, "\\1", desc);
    vers <- gsub(pattern, "\\2", desc);
    vers <- gsub("[()]", "", vers);
    names(vers) <- pkgs;
    dups <- duplicated(pkgs);
    vers <- vers[!dups];
    version <- vers["BiocParallel"];
  }
  R.utils::use("BiocParallel", version=version, ...);
} # .useBiocParallel()


############################################################################
# HISTORY:
# 2015-09-05
# o Add support for dsApply(..., .parallel="future").
# 2015-01-05
# o CLEANUP: Using requireNamespace() instead of require() internally.
# 2014-08-07
# o BUG FIX: dsApply() for GenericDataFileSet would coerce argument
#   'verbose' to logical before applying the function.
# 2014-04-19
# o dsApply(..., .parallel="none") would lower the verbose threshold
#   before applying the function.
# 2014-03-30
# o Added dsApplyInPairs().
# o Now dsApply(..., .parallel="BatchJobs") also returns the results.
# o Moved to the R.filesets package (from aroma.seq).
# 2014-01-24
# o Now argument 'IDXS' can also be an index vector, which is then
#   treated as as.list(IDXS).
# 2014-01-04
# o CONSISTENCY: Now dsApply() returns a named vector with names
#   corresponding to the names of the processed file items.
# 2013-11-22
# o Added argument 'IDXS' to dsApply() allowing to process not only
#   individuals files but also lists of individuals files.
# 2013-11-02
# o Adding support for distributed processing via 'BiocParallel'.
# 2013-11-01
# o Made the code more generic.
# 2013-09-28
# o ROBUSTNESS: Now .getBatchJobRegistryId() adds process ID ("pid") and
#   a timestamp to the registry path to make it even more unique.
# o BUG FIX: dsApply() for GenericDataFileSet when executed via the
#   'BatchJobs' methods would not allow commas in the work directory
#   among other directories.
# 2013-08-31
# o Now .useBatchJobs() utilizes R.utils::use().
# o Added dsApply() for GenericDataFileSet.  Added Rdocs.
# 2013-08-26
# o Added .getBatchJobRegistry() and .getBatchJobRegistryId().
# o Added .usePackage() and .useBatchJobs().
# o Created.
############################################################################
