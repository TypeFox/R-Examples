
compareTranscriptAndOutput <- function(name, tests, results, verbose=TRUE) {
    x <- list()
    # return a RtTestSetResults objects (a list of lists), each having components:
    #   status: one of "error", "warning", "info", "ok"
    #   msg: an informative message (multi-line as char vector)
    is.empty.command <- function(case) {
        return(is.element("comment", names(case)) && identical(case$comment, "> "))
    }
    ## see if we can strip empty commands from the end of the longer set
    while (length(results) > length(tests) && is.empty.command(results[[length(results)]]))
        results <- results[-length(results)]
    while (length(tests) > length(results) && is.empty.command(tests[[length(tests)]]))
        tests <- tests[-length(tests)]
    if (length(tests) != length(results)) {
        msg <- "Number of commands in tests and results files differs!"
        haltedEarly <- FALSE
        if (length(tests) > length(results)) {
            msg <- c(msg, "  tests file has more commands, first extra one is:",
                     paste("  ", tests[[length(results)+1]]$input, sep=""))
            xtraTests <- tests[-seq(along=results)]
            xtraResults <- list()
            tests <- tests[seq(along=results)]
            if (any(regexpr("execution halted", results[[length(results)]]$output, ignore.case=TRUE)>0)) {
                haltedEarly <- TRUE
                msg <- c(msg, "Execution of tests was halted early")
            }
        } else {
            msg <- c(msg, "  results file has more commands, first extra one is:",
                     paste("  ", results[[length(tests)+1]]$input, sep=""))
            xtraResults <- results[-seq(along=tests)]
            xtraTests <- list()
            results <- results[seq(along=tests)]
        }
        if (verbose)
            cat(paste("Error in", name), msg, sep="\n")
        # need double nesting of list because of unlist( , recursive=FALSE) below
        x <- c(x, list(list(list(msg=msg, status="error"))))
    }
    x <- c(x, mapply(tests, results, seq(along=tests), FUN=function(test, result, testnum) {
        ## check that input matches
        res <- list()
        if (!identical(all.equal(test$input, result$input), TRUE)) {
            res <- list(list(status="warning", testnum=testnum,
                             msg=c("Internal inconsistency: mismatch between test input & input from transcript",
                                   "Test input:", test$input,
                                   "Transcript input:", result$input)))
        }
        res <- c(res, list(compareSingleTest(test$input, test$control, test$output, test$comment, test$garbage, result$output, testnum, name, verbose)))
        return(res)
    }, SIMPLIFY=FALSE, USE.NAMES=FALSE))
    x <- unlist(x, use.names=FALSE, recursive=FALSE)
    if (verbose)
        cat("\n")
    class(x) <- "RtTestSetResults"
    attr(x, "testname") <- name
    x
}
