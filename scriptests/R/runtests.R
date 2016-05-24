runtests <- function(pkg.dir=getOption("scriptests.pkg.dir", "pkg"),
                     pattern=".*", file=NULL,
                     full=FALSE, dir=TRUE,
                     clobber=FALSE, output.suffix=NULL, console=FALSE,
                     ...,
                     verbose=TRUE, envir=globalenv(), subst=NULL,
                     path=getOption("scriptests.pkg.path", default=getwd())) {
    if (!interactive() && basename(getwd())=="tests") {
        if (nargs() != 0)
            warning("runtests() is for interactive use - use runScripTests() in tests/runtests.R - calling runScripTests() and ignoring all arguments supplied to runtests()")
        # Looks like we're being called by R CMD check (because runtests()
        # instead of runScripTests() was put in tests/runtests.R.
        status <- runScripTests()
        return(status)
    }
    supplied.path <- !missing(path)
    supplied.pkg.dir <- !missing(pkg.dir)
    orig.path <- path

    if (regexpr("^(/|\\\\|[a-zA-Z]:)", pkg.dir) > 0) {
        pkg.dir.path <- pkg.dir
        pkg.dir <- basename(pkg.dir)
        path.used <- dirname(pkg.dir)
    } else {
        for (path.used in path.expand(strsplit(path, ';')[[1]])) {
            pkg.dir.path <- pkg.path(path.used, pkg.dir)
            if (file.exists(pkg.dir.path))
                break
        }
    }
    if (!file.exists(pkg.dir.path))
        stop("cannot find package directory ", pkg.dir.path, " using path='", path.used, "'")
    pkg.name <- read.pkg.name(pkg.dir.path, pkg.dir)
    if (!is.null(file) && pattern!=".*")
        stop("cannot supply both 'file' and 'pattern'")
    if (is.null(subst) && !full) {
        # remove text "package:::" if "package" is not attached
        if (!is.element(paste("package:", pkg.name, sep=""), search())) {
            subst <- rep("", 2)
            names(subst) <- c(paste("\\b", pkg.name, ":::?", sep=""), paste("\\blibrary\\( *['\"]?", pkg.name, "['\"]? *\\)", sep=""))
            cat("* Package '", pkg.name, "' is not loaded as a package; will remove ", paste('"', names(subst), '"', collapse=", ", sep=""), " in tests\n", sep="")
        }
    } else if (is.logical(subst) && !subst) {
        # Don't do any substitution
        subst <- NULL
    } else if (!is.character(subst) && !is.null(subst)) {
        stop("subst must be NULL, F, or a character vector")
    }
    if (full) {
        if (!is.null(output.suffix))
            stop("can only supply output.suffix= argument when full==FALSE")
        if (!is.null(file))
            stop("can only supply file= argument when full==FALSE")
    }
    # Prepare to change directory if needed
    cwd <- getwd()
    test.dir <- file.path(pkg.path(path.used, pkg.dir), "tests")
    # convert to absolute path because this function changes working directories,
    test.dir <- normalizePath(test.dir)
    if (!file.exists(test.dir))
        stop("test source directory ", test.dir, " does not exist")
    dir.isStandard <- FALSE
    if (is.logical(dir)) {
        if (dir) {
            if (full)
                dir <- paste(pkg.name, ".Rcheck/tests", sep="")
            else
                dir <- paste(pkg.dir, ".tests", sep="")
            dir.isStandard <- TRUE
        } else {
            dir <- NULL
        }
    }
    # Need this later, relative to the current path
    if (full) {
        check.dirs <- paste(c(pkg.name, pkg.dir), ".Rcheck", sep="")
        if (!any(file.exists(check.dirs)))
            stop("expected to find installed library ", pkg.name, " in ",
                 paste("'", check.dirs, "'", sep="", collapse=" or "),
                 " but neither of those directories exists")
        check.dirs <- check.dirs[file.exists(check.dirs)]
        if (!any(file.exists(file.path(check.dirs, pkg.name))))
            stop("expected to find installed library in ", file.path(check.dirs, pkg.name), " but that directory doesn't exist")
        check.dirs <- check.dirs[file.exists(file.path(check.dirs, pkg.name))]
        # if there's more than one, choose the one with the most recent modification time
        mt <- file.info(file.path(check.dirs, pkg.name))[,"mtime"]
        if (length(check.dirs) > 1)
            check.dirs <- check.dirs[which.max(mt)]
        cat("* Using package in '", file.path(check.dirs, pkg.name), "' for running tests\n", sep="")
        # This code relies on normalizePath converting to an absolute path
        r.libs.site.orig <- Sys.getenv("R_LIBS_SITE")[[1]]
        r.libs.site <- unique(sapply(c(strsplit(r.libs.site.orig, split=";")[[1]], check.dirs), normalizePath))
        Sys.setenv("R_LIBS_SITE"=paste(r.libs.site, collapse=";"))
        on.exit(Sys.setenv("R_LIBS_SITE"=r.libs.site.orig), add=TRUE)
    }
    if (path.used=='.')
        path.used <- getwd()
    if (!is.null(dir)) {
        if (file.exists(dir)) {
            if (!clobber && !dir.isStandard)
                stop("dir for running tests '", dir, "' already exists and is non-standard - supply clobber=TRUE to overwrite")
            if (verbose)
                cat("* Removing old tests directory ", dir, "\n", sep="")
            unlink(dir, recursive=TRUE)
        }
        if (!dir.create(dir, recursive=TRUE))
            stop("failed to create directory: ", dir)
        if (verbose)
            cat("* Copying ", test.dir, " to ", dir, "\n", sep="")
        for (f in list.files(test.dir))
            file.copy(file.path(test.dir, f), dir, recursive=TRUE)
        existing.files <- list.files()
        if (verbose)
            cat("* Setting working directory to ", dir, "\n", sep="")
        setwd(dir)
        on.exit(setwd(cwd), add=TRUE)
    }
    if (!full) {
        if (!is.null(output.suffix) && length(output.suffix)!=1) {
            ## Functions installed in on.exit() do not always seem to reliably
            ## run after an error in a called function, so make habit of calling
            ## them manually before an error.
            eval(sys.on.exit())
            stop("length(output.suffix)!=1")
        }
        if (length(list(...)))
            warning("ignoring extra arguments when full=FALSE: ", paste(names(list(...)), collapse=", "))
    }

    if (supplied.pkg.dir)
        options("scriptests.pkg.dir"=pkg.dir)
    if (supplied.path)
        options("scriptests.pkg.path"=orig.path)

    if (!full) {
        old.options <- options(width=80, warn=1)
        on.exit(options(old.options), add=TRUE)
        res <- runTestsHereFast(pattern=pattern, pkg.dir=pkg.dir, pkg.name=pkg.name, file=file, verbose=verbose, envir=envir, subst=subst, path=path.used)
        attr(res, "dir") <- dirname(names(res)[1])
        names(res) <- basename(names(res))
        if (console || (!is.null(output.suffix) && !(is.logical(output.suffix) && !output.suffix)))
            dumprout(res, output.suffix=output.suffix, console=console, clobber=clobber)
        return(invisible(res))
    } else {
        status <- runScripTests(..., quit=FALSE, subst=subst, pattern=pattern)
        new.files <- setdiff(list.files(), existing.files)
        if (FALSE && length(new.files)) {
            cat("* Removing ", length(new.files), " new files: ", paste(new.files, collapse=", "), "\n", sep="")
            file.remove(new.files)
        }
        return(status)
    }
}
