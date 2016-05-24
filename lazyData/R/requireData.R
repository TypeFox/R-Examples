requireData <- function(package = stop("you must specify a package"),
                        lib.loc = NULL, quietly = TRUE, character.only = FALSE,
                        warn.conflicts = TRUE, reallyQuietly = TRUE, ...) {
  if(!character.only) {
    pkg <- substitute(package)
    if(!is.character(pkg)) pkg <- deparse(pkg)
  } else pkg <- as.character(package)

  stopifnot(length(pkg) == 1)

  ## check if the package is on, or can be put on, the search path:
  s0 <- search()
  oldWarn <- options(warn = -1)
  on.exit(options(oldWarn))
  OK <- if(reallyQuietly) {
    suppressPackageStartupMessages(require(package = pkg, lib.loc = lib.loc,
                                           quietly = TRUE, warn.conflicts = FALSE,
                                           character.only = TRUE))
  } else {
    require(package = pkg, lib.loc = lib.loc, quietly = quietly,
            warn.conflicts = warn.conflicts, character.only = TRUE)
  }
  options(oldWarn)
  if(!OK) {
    warning("no valid package called ", sQuote(pkg), " found")
    return(invisible(FALSE))  ### package not installed
  }

  ## take care of any extra packages as well
  packages <- unique(c(pkg, sub("^package:", "",
                            grep("^package:", setdiff(search(), s0), value = TRUE))))

  for(package in packages)
      ## check if it needs to have lazy data
      if(file.exists(f <- system.file("data", package = package)) &&
         file.info(f)$isdir &&
         !file.exists(file.path(f, "Rdata.rds"))) {

        ## check if the package has any data sets
        d <- utils::data(package = package)$results
        if(nrow(d) == 0) next  ## no data sets to expose

        d <- d[, "Item"]  ## names of data sets in form 'objName (dataName)'
        objName <- sub(" .*$", "", d)
        datName <- sub("^.*\\(", "", sub("\\)$", "", d))

        ## remove any prior attachments of the datasets:
        searchDataName <- paste0("datasets:", package)
        while(searchDataName %in% search()) detach(searchDataName,
                                                   character.only = TRUE)

        ## find out where the package now sits on the search path:
        searchPkgName <- paste0("package:", package)
        pkgIndex <- match(searchPkgName, search(), nomatch = 1)

        ## find out where to put the datasets and check OK:
        pos <- pkgIndex + 1

        ## set up the datasets entry and populate it with 'promises':
        env <- attach(NULL, pos = pos, name = searchDataName)
        attr(env, "path") <- attr(as.environment(pkgIndex), "path")

        for(i in seq_along(d))
            eval(substitute(delayedAssign(OBJ, {
              data(DAT, package = PKG,
                   envir = as.environment(match(SEARCHDATANAME, search())))
              get(OBJ)
            }, eval.env = env, assign.env = env),
                            list(OBJ = objName[i],
                                 SEARCHDATANAME = searchDataName,
                                 PKG = package,
                                 DAT = datName[i])))
      }
  invisible(TRUE)
}
