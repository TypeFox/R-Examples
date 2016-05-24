# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sourceTestScript <- function(pathname, devel=FALSE, ...) {
  # Argument 'devel':
  devel <- Arguments$getLogical(devel);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Enable support for sibling root paths?
  if (devel) {
    key <- "devel/dropRootPathTags";
    oopts <- setOption(aromaSettings, key, TRUE);
    on.exit({
      setOption(aromaSettings, key, oopts);
    }, add=TRUE);
    # Print settings
    print(getOption(aromaSettings, key));
  }

  # Use special file cache for testing
  library("R.cache");
  opath <- getCacheRootPath();
  on.exit({
    setCacheRootPath(opath);
  }, add=TRUE);
  setCacheRootPath(".Rcache");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cat("** PATHNAME: ", pathname, "\n", sep="");
  tryCatch({
    pathname <- Arguments$getReadablePathname(pathname);
    source(pathname, echo=TRUE);
    cat("** PATHNAME DONE: ", pathname, "\n", sep="");
  }, error = function(ex) {
    cat("************************************************\n");
    cat("** ", rep(c(" ER", "ROR "), times=6), " **\n", sep="");
    print(ex);
    print(sessionInfo());
    cat("************************************************\n");
    cat("** PATHNAME FAILED: ", pathname, "\n", sep="");
  });
  gc();
} # sourceTestScript()


sourceTestSet <- function(path, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each test script...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnames <- list.files(path=path, pattern="[.]R$", full.names=TRUE);
  cat("TEST SCRIPTS: \n");
  print(pathnames);

  for (pathname in pathnames) {
    if (regexpr("hetero", pathname) != -1)
      next;
    if (regexpr("expectile", pathname) != -1)
      next;
    sourceTestScript(pathname, ...);
  } # for (pathname ...);
} # sourceTestSet()


sourceTestGroups <- function(groups, pattern=".*", order=c("auto", "asis", "reverse", "random"), nbrOfSets=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'groups':
  groups <- Arguments$getCharacters(groups);

  # Argument 'pattern':
  pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'order':
  order <- match.arg(order);

  # Argument 'nbrOfSets':
  if (!is.null(nbrOfSets)) {
    nbrOfSets <- Arguments$getInteger(nbrOfSets, range=c(1,Inf));
    if (order == "auto") order <- "random";
  }

  if (order == "auto") {
    order <- "asis";
  }

  rootPath <- Arguments$getReadablePath("testScripts/");
  for (group in groups) {
    path <- file.path(rootPath, group);
    path <- Arguments$getReadablePath(path);

    # Find all test sets
    pathnames <- list.files(path=path, recursive=TRUE, full.names=TRUE);
    paths <- dirname(pathnames);
    paths <- unique(paths);
    paths <- paths[sapply(paths, FUN=isDirectory)];

    # Filter
    paths <- grep(pattern, paths, value=TRUE);
    cat("Filtered test sets:\n");
    print(paths);

    if (order == "reverse") {
      paths <- rev(paths);
    } else if (order == "random") {
      paths <- sample(paths);
    }
    cat("All test sets:\n");
    print(paths);

    if (!is.null(nbrOfSets)) {
      cat("Subsetted test sets:\n");
      paths <- paths[1:nbrOfSets];
      print(paths);
    }

    # For each test set...
    for (path in paths) {
      sourceTestSet(path, ...);
    }
  } # for (groups ...);
} # sourceTestGroups()


launchTestGroups <- function(groups="system", ...) {
  sourceTestGroups(groups=groups, ...);
} # launchTestGroups()



############################################################################
# HISTORY:
# 2012-10-19
# o BUG FIX: The temporary cache directory was set to 'TRUE'.
# 2012-09-14
# o Created.
############################################################################
