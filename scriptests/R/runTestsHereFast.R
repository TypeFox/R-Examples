runTestsHereFast <- function(pattern=".*",
                             pkg.dir=getOption("scriptests.pkg.dir"),
                             pkg.name=NULL,
                             file=NULL,
                             verbose=TRUE, envir=globalenv(), subst=NULL,
                             test.suffix=".Rt",
                             path=getOption("scriptests.pkg.path", default=getwd())) {
    # This does the similar work as runScripTests()/.runPackageTests(),
    # with these differences:
    #
    #   (1) tests are run in the current directory rather than creating
    #       a copy of the package 'tests' directory and doing setwd() on it
    #   (2) all test code is run in this R session (runScripTests() runs
    #       each file in a different R session)
    #   (3) doesn't read the CONFIG file
    #   (4) use of scriptests initialize/diff/finalize is hardwired in here
    #   (5) output is captured using evalCapture() instead of reading it
    #       from a transcript
    pkg.dir.path <- pkg.path(path, pkg.dir)
    if (is.null(pkg.name))
        pkg.name <- read.pkg.name(pkg.dir.path, pkg.dir)
    if (!is.null(file)) {
        files <- file
        if (!all(i <- file.exists(files))) {
            warning("ignoring non-existant files ", paste(files[!i], collapse=", "))
            files <- files[i]
        }
    } else {
        if (nchar(test.suffix))
            test.suffix <- gsub("^\\.", "\\\\.", test.suffix)
        if (regexpr("\\$$", pattern) < 1
            && regexpr(paste(test.suffix, "\\$?$", sep=""), pattern, ignore.case=T) < 1) {
            pattern <- paste(pattern, ".*", test.suffix, sep="")
            if (regexpr("\\$$", test.suffix) < 1)
                pattern <- paste(pattern, "$", sep="")

        }
        files <- list.files('.', pattern=pattern, full.names=TRUE, ignore.case=TRUE)
        if (length(files)==0)
            stop("no files matched the pattern '", pattern, "' in ", getwd())
    }
    allres <- list()
    for (file in files) {
        if (verbose)
            cat("* Running tests in", file)
        tests <- parseTranscriptFile(file, subst=subst)
        if (verbose)
            cat(" (read", length(tests), "chunks)\n")
        res <- lapply(seq(along=tests), function(i) {
            test <- tests[[i]]
            if (is(test$expr, "try-error"))
                actual <- as.character(test$expr)
            else
                actual <- evalCapture(test$expr, envir)
            res <- compareSingleTest(test$input, test$control, test$output, test$comment, test$garbage, actual,
                                     i, file, verbose=verbose)
            res$comment <- test$comment
            res$transcript <- c(test$input, test$control, actual)
            res$target <- c(test$output)
            res
        })
        class(res) <- "RtTestSetResults"
        attr(res, "testname") <- file
        if (verbose) {
            cat("\n")
            print(summary(res))
        }
        allres[[file]] <- res
    }
    class(allres) <- "RtTestSetResultsList"
    attr(allres, "dir") <- getwd() # file.path(pkg.path(path, pkg.dir), "tests")
    attr(allres, "pattern") <- pattern
    if (length(allres)>1)
        print(summary(allres))
    invisible(allres)
}

