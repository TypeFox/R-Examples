print.ctv <- function(x, packagelist = TRUE, ...)
{
  cat(paste("\nCRAN Task View\n--------------\nName:       ",
            x$name, "\nTopic:      ",
	    x$topic, "\nMaintainer: ",
	    x$maintainer,
	    if(!is.null(x$email)) sprintf("\nContact:    %s", x$email),
	    sprintf("\nVersion:    %s", x$version),
	    "\n",
	    ifelse(is.null(x$repository), "", paste("Repository: ", x$repository, sep = "")),
	    "\n", sep = ""))
  if(packagelist) {
    pkgs <- paste(strwrap(paste(x$packagelist$name, ifelse(x$packagelist$core, "*", ""), sep = "", collapse = ", "),
      width = getOption("width"), prefix = "            "), collapse = "\n")
    substr(pkgs, 1, 9) <- "Packages:"
    cat(pkgs)
    cat(ifelse(any(x$packagelist$core), "\n            (* = core package)\n\n", "\n"))
  }
  invisible(x)
}

print.ctvlist <- function(x, packagelist = FALSE, ...)
{
  cat("\nCRAN Task Views\n")
  n <- length(x)
  if(n < 1)
    cat("no views found\n")
  else {
    for(i in seq(along = x)) {    
      cat(paste("---------------\nName: ",
                x[[i]]$name, "\nTopic: ",
	        x[[i]]$topic, "\nMaintainer: ",
	        x[[i]]$maintainer, "\n",
	        ifelse(is.null(x[[i]]$repository), "", paste("Repository: ", x[[i]]$repository, "\n", sep = "")),
	        sep = ""))
      if(packagelist) {
        cat(paste("Packages:",
            paste(x[[i]]$packagelist$name, ifelse(x[[i]]$packagelist$core, "*", ""), sep = "", collapse = ", "),
            "\n"))
        cat(ifelse(any(x[[i]]$packagelist$core), "(* = core package)\n", ""))
      }
    }
    cat("\n")
  }
  invisible(x)
}

available.views <- CRAN.views <- function(repos = NULL, ...)
{
  ## getOption("repos") replaces getOption("CRAN") from 2.1.0 on
  if(is.null(repos)) repos <- ifelse(is.null(getOption("repos")), getOption("CRAN"), getOption("repos"))

  if("@CRAN@" %in% repos && interactive()) {
      cat(gettext("--- Please select a CRAN mirror for use in this session ---\n"))
      flush.console()
      chooseCRANmirror()
      m <- match("@CRAN@", repos)
      nm <- names(repos)
      repos[m] <- getOption("repos")["CRAN"]
      if(is.null(nm)) nm <- rep("", length(repos))
      nm[m] <- "CRAN"
      names(repos) <- nm
  }
  if("@CRAN@" %in% repos) stop("trying to use CRAN without setting a mirror")

  contriburl <- paste(repos, "/src/contrib", sep = "")
  rval <- list()
  
  for(i in seq(along.with = contriburl)) {
    ## load Views.rds from repository
    viewurl <- gzcon(url(paste(contriburl[i], "Views.rds", sep = "/"), open = "rb"))
    x <- suppressWarnings(try(readRDS(viewurl), silent = TRUE))
    if(inherits(x, "try-error")) next else close(viewurl)

    ## add repository information    
    for(j in seq(along = x)) x[[j]]$repository <- repos[i]
    rval <- c(rval, x)
  }
  class(rval) <- "ctvlist"
  return(rval)
}

.get_pkgs_from_ctv_or_repos <- function(views, coreOnly = FALSE, repos = NULL)
{
  ## views is already a "ctv"
  if(inherits(views, "ctv")) {
    views <- list(views)
    class(views) <- "ctvlist"
  }
  
  ## views is a string that needs to be turned to "ctv"
  if(!inherits(views, "ctvlist")) {
    ## get CRAN views and extract names of available views
    cranviews <- available.views(repos = repos)
    availnames <- sapply(seq(along = cranviews), function(i) cranviews[[i]]$name)

    whichviews <- lapply(views, function(z) {
      rval <- which(z == availnames)
      if(length(rval) > 0) rval[1] else numeric(0)
    })
    unavail <- which(sapply(whichviews, length) < 1)
    if(length(unavail) > 0) warning(paste("CRAN task view", views[unavail], "not available", collapse = "\n"))
    views <- cranviews[as.vector(unlist(whichviews))]
    class(views) <- "ctvlist"
  }

  ## compute list of unique packages for each repository
  coreOnly <- rep(coreOnly, length.out = length(views)) 
  pkgs <- lapply(seq(along = views), function(i) {
    pkgs <- if(!coreOnly[i]) views[[i]]$packagelist$name
      else views[[i]]$packagelist$name[views[[i]]$packagelist$core]
    cbind(pkgs, views[[i]]$repository)
  })
  pkgs <- do.call("rbind", pkgs)
  pkgs <- tapply(pkgs[,1], pkgs[,2], function(x) sort(unique(x)))
  
  return(pkgs)
}

install.views <- function(views, coreOnly = FALSE, repos = NULL, ...)
{
  pkgs <- .get_pkgs_from_ctv_or_repos(views = views, coreOnly = coreOnly, repos = repos)
  for(i in seq(along = pkgs)) install.packages(pkgs[[i]], repos = names(pkgs)[i], ...)    
  invisible()
}

download.views <- function(views, destdir, coreOnly = FALSE, repos = NULL, ...)
{
  pkgs <- .get_pkgs_from_ctv_or_repos(views = views, coreOnly = coreOnly, repos = repos)
  for(i in seq(along = pkgs)) download.packages(pkgs[[i]], destdir = destdir, repos = names(pkgs)[i], ...)    
  invisible()
}

update.views <- function(views, coreOnly = FALSE, repos = NULL, lib.loc = NULL, filters = NULL, ...)
{
  pkg_list <- .get_pkgs_from_ctv_or_repos(views = views, coreOnly = coreOnly, repos = repos)
    
  for(i in seq(along = pkg_list)) {
    ## package and repos name(s)
    pkgs <- pkg_list[[i]]
    repos <- names(pkg_list)[i]
        
    ## query available packages
    apkgs <- available.packages(contriburl = contrib.url(repos), filters = filters)

    ## compute intersection
    unavail <- which(!(pkgs %in% as.character(apkgs[, 1L])))
    if(length(unavail) > 0L) {
      warning(sprintf("The following packages are not available: %s", paste(pkgs[unavail], collapse = ", ")))
      pkgs <- pkgs[-unavail]
    }
    apkgs <- apkgs[pkgs, , drop = FALSE]

    ## query installed packages
    ipkgs <- installed.packages(lib.loc = lib.loc)
    ipkgs <- ipkgs[which(ipkgs[, "Package"] %in% pkgs), , drop = FALSE]
  
    ## determine which packages need to be updated
    if(NROW(ipkgs) > 0L) {    
      not_installed <- which(!(pkgs %in% ipkgs[, 1L]))
    
      cpkgs <- if(length(not_installed) > 0L) pkgs[-not_installed] else pkgs
      get_highest_version <- function(x)
        as.character(max(package_version(ipkgs[which(ipkgs[, 1L] %in% x),"Version"])))
      not_uptodate <- which(package_version(apkgs[cpkgs, "Version"]) >
        package_version(sapply(cpkgs, get_highest_version)))
      pkgs <- sort(c(pkgs[not_installed], cpkgs[not_uptodate]))
    }

    ## install packages required
    apkgs <- apkgs[pkgs, , drop = FALSE]
    if(NROW(apkgs) > 0L) install.packages(apkgs[, 1L], repos = repos, lib = lib.loc, ...)
  }
  
  invisible()
}
