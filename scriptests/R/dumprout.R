dumprout <- function(res = .Last.value, output.suffix = ".Rout.tmp", verbose = TRUE, console = FALSE, files = !console, clobber = identical(output.suffix, ".Rout.tmp"), level=c("error", "all", "info", "warning")) {
    level <- match.arg(level)
    if (inherits(res, "RtTestSetResults")) {
        # make this single test result look like a list of results
        res <- list(res)
        names(res) <- basename(attr(res[[1]], "testname"))
        class(res) <- "RtTestSetResultsList"
    }
    if (!inherits(res, "RtTestSetResultsList"))
        stop("supplied argument is not a return value of runtests()")
    if (length(output.suffix)!=1)
        stop("length(output.suffix)!=1")
    if (is.character(output.suffix) && output.suffix!="") {
        if (substring(output.suffix, 1, 1)!=".")
            output.suffix <- paste(".", output.suffix, sep="")
    } else {
        output.suffix <- ".Rout.tmp"
    }
    for (i in seq(along=res)) {
        outfile <- basename(attr(res[[i]], "testname"))
        if (level != "all") {
            cc <- summary(res[[i]])$counts
            if (   (level == "error" && cc["error"] == 0)
                || (level == "warning" && all(cc[c("warning", "error")] == 0))
                || (level == "info" && all(cc[c("info", "warning", "error")] == 0))) {
                if (verbose)
                    cat("* Skipping ", outfile, ": no notifications at level '", level, "'",
                        if (level!="error") " or above", "\n", sep="")
                next
            }
        }
        if (!files) {
            print(res[[i]], transcript=TRUE)
        } else {
            if (outfile != "") {
                outfile <- gsub("\\.Rt", output.suffix, outfile)
                if (file.exists(outfile) && !clobber) {
                    warning("file ", outfile, " exists already, specify clobber=TRUE to overwrite")
                } else {
                    if (verbose)
                        cat("* Writing transcript of actual output to ", outfile, "\n", sep="")
                    sink(outfile)
                    print(res[[i]], transcript=TRUE)
                    sink()
                }
            } else {
                warning("no filename for output from test file ", i, ": nowhere to write it")
            }
        }
    }
    invisible(res)
}
