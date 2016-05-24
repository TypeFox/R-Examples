# runScripTests() is intended to be called from a .R file in the tests directory
# e.g., like this:
#    library(ScripTests)
#    runScripTests()

runScripTests <- function(..., initializeFun = Quote(initializeTests()),
                            finalizeFun = Quote(summarizeTests()),
                            diffFun = Quote(ScripDiff()), subst=NULL,
                            pattern=NULL, quit=TRUE)
{
    # When run by "R CMD check", no output directly from the tests will appear on the
    # console.  The only output the user sees will be the last 13 lines from the first
    # .fail file that "R CMD check" files after getting a non-zero status result from
    # running a .R file.
    # Try to work out whether we were called from a .R file, and if so what was its name,
    # so that we can make the name of the file with the transcript of detailed comparisons
    # appear at the end of the .fail file so that it will be displayed to the user.
    test.transcript.file <- character(0)
    run.from <- character(0)
    if (!is.na(i <- match("-f", commandArgs()))) {
        run.from <- basename(commandArgs()[i+1])
        test.transcript.file <- paste(commandArgs()[i+1], "out", sep="")
        test.transcript.file <- paste(c(rev(rev(strsplit(gsub("\\\\", "/", getwd(), perl=TRUE), "/")[[1]])[1:2]), test.transcript.file), collapse=.Platform$file.sep)
    }
    cat("\n");
    status <- .runPackageTests(..., initializeFun=initializeFun, finalizeFun=finalizeFun,
                               diffFun=diffFun, run.preexisting.R.files=FALSE, subst=subst,
                               pattern=pattern, run.from=run.from)
    if (length(test.transcript.file))
        cat("\nSee ", test.transcript.file, (if (status) ".fail"), " for a transcript of test comparisons\n", sep="")
    fail.files <- list.files(pattern="\\.Rout\\.fail$")
    fail.files <- setdiff(fail.files, basename(test.transcript.file))
    if (length(fail.files))
        cat("Look for clues in ", if (length(fail.files)==1) "this file" else "these files",
            if (length(test.transcript.file)) " too: ",
            paste("'", fail.files, "'", collapse=", ", sep=""), "\n", sep="")
    if (status && !(Sys.getenv("SCRIPTESTS13") %in% c("0", "F", "FALSE"))
        && exists("test-summary.fail", envir=test.status.env, inherits=FALSE)) {
        # output 13 lines of low-level error comparison for R CMD check to output
        testResults <- get("test-summary.fail", envir=test.status.env)
        if (length(testResults)>1) {
            oneBadFile <- gsub(":$", ".log", testResults$name[length(testResults$name)-1])
            if (file.exists(oneBadFile)) {
                lines <- readLines(oneBadFile)
                if (length(lines)) {
                    cat("***** 13 lines of low-level output here b/c R CMD check wants to display just 13 lines, turn this off by setting environment variable SCRIPTESTS13=0\n")
                    lines <- lines[seq(max(1, length(lines)-12), length(lines))]
                    cat(lines, sep="\n")
                    if (length(test.transcript.file))
                        cat("See ", test.transcript.file, ".fail for more info\n", sep="")
                }
            }
        }
    }
    if (quit && !interactive())
        q("no", status = status)
    invisible(NULL)
}
