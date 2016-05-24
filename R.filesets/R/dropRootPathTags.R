setMethodS3("dropRootPathTags", "default", function(path, depth, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getCharacter(path);

  # Argument 'depth':
  depth <- Arguments$getInteger(depth, range=c(0,Inf));

  rootPath <- getParent(path, depth=depth);
  if (length(rootPath) == 0) {
    throw("Argument 'depth' (", depth, "') is too large for this path: ", path);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Dropping tags from root path");

  verbose && cat(verbose, "Path: ", path);

  verbose && cat(verbose, "Root path: ", rootPath);

  rootRootPath <- dirname(rootPath);
  rootPath <- basename(rootPath);
  rootPath <- gsub(",.*", "", rootPath);
  if (rootRootPath != ".") {
    rootPath <- file.path(rootRootPath, rootPath);
  }
  verbose && cat(verbose, "Root path without tags: ", rootPath);

  subdirs <- sapply(seq_len(depth), FUN=function(d) {
    basename(getParent(path, depth=d-1L));
  });
  subdirs <- rev(subdirs);
  subdirs <- do.call(file.path, args=as.list(subdirs));
  verbose && cat(verbose, "Subdirectories: ", subdirs);

  path <- file.path(rootPath, subdirs);
  verbose && cat(verbose, "Path without root-path tags: ", path);

  verbose && exit(verbose);

  path;
}, protected=TRUE) # dropRootPathTags()

############################################################################
# HISTORY:
# 2011-02-24
# o Added dropRootPathTags().
############################################################################  
