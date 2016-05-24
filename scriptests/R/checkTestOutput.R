checkTestOutput <- function(rtIn, rtSave=paste(rtIn, ".save", sep=""), debug=TRUE) {
    rtOut <- gsub("\\.Rt$", ".Rout", rtIn, perl=TRUE)
    rtFail <- gsub("\\.Rt$", ".Rt.fail", rtIn, perl=TRUE)
    rtLog <- paste(rtIn, ".log", sep="")
    if (!file.exists(rtIn)) {
        msg <- paste("checkTestOutput: cannot find original test file '", rtIn, "' in '", getwd(), "'\n", sep="")
        cat(file=rtFail, msg)
        cat(file=stdout(), msg)
        return(NULL)
    }
    if (!file.exists(rtSave)) {
        msg <- paste("checkTestOutput: cannot find saved-test-object file '", rtSave, "' in '", getwd(), "'\n", sep="")
        cat(file=rtFail, msg)
        cat(file=stdout(), msg)
        return(NULL)
    }
    if (!file.exists(rtOut)) {
        msg <- paste("checkTestOutput: cannot actual test output file '", rtOut, "' in '", getwd(), "'\n", sep="")
        cat(file=rtFail, msg)
        cat(file=stdout(), msg)
        return(NULL)
    }
    # sink(file=stdout())
    if (debug)
        cat("  * Loading saved transcript object from file \"", rtSave, "\" ...\n", sep="", file=stdout())

    testObjName <- load(file=rtSave, envir=as.environment(-1))
    if (testObjName[1] != "tests")
        tests <- get(testObjName[1])
    cat("  * Parsing actual test output from file \"", rtOut, "\" ...\n", sep="", file=stdout())
    resList <- parseTranscriptFile(rtOut, ignoreUpToRegExpr="> # End of scriptests preamble", ignoreAfterRegExpr="> # End of RtTests output")
    res <- compareTranscriptAndOutput(sub(".Rout", ".Rt", rtOut), tests, resList, verbose=TRUE)
    res.summary <- summary(res)
    print(res.summary)
    # sink()
    sink(rtLog)
    print(res)
    print(res.summary)
    sink()
    # Incrementally add to the summary in test-summary.txt, updating the total
    testResultsFile <- "test-summary.txt"
    # Test summaries are lines like this:
    #   1        2 3     4    5 6       7 8        9   10 11
    #   plus.Rt: 4 tests with 0 errors, 0 warnings and 0 messages
    if (file.exists(testResultsFile)) {
        testResults <- scan("test-summary.txt", quiet=T, what=list(name="", ntests=0, NULL, NULL, nerr=0, NULL, nwarn=0, NULL, NULL, nmess=0, NULL), fill=T)
        if (length(testResults$name)>0 && testResults$name[length(testResults$name)]=="total:")
            testResults <- lapply(testResults, "[", -length(testResults$name))
    } else {
        testResults <- list(name=character(0), ntests=numeric(0), c1=NULL,
                                  c2=NULL, nerr=numeric(0), c3=NULL,
                                  nwarn=numeric(0), c4=NULL, c5=NULL,
                                  nmess=numeric(0), c6=NULL)
    }
    i <- match(paste(rtIn, ":", sep=""), testResults$name, nomatch=length(testResults$name)+1)
    testResults$name[i] <- paste(rtIn, ":", sep="")
    testResults$ntests[i] <- res.summary$n
    testResults$nerr[i] <- res.summary$counts["error"]
    testResults$nwarn[i] <- res.summary$counts["warning"]
    testResults$nmess[i] <- res.summary$counts["info"]
    # add a line for the totals
    i <- length(testResults$name)+1
    testResults$name[i] <- "total:"
    testResults$ntests[i] <- sum(testResults$ntests, na.rm=T)
    testResults[[3]] <- "tests"
    testResults[[4]] <- "with"
    testResults$nerr[i] <- sum(testResults$nerr, na.rm=T)
    testResults[[6]] <- "errors,"
    testResults$nwarn[i] <- sum(testResults$nwarn, na.rm=T)
    testResults[[8]] <- "warnings"
    testResults[[9]] <- "and"
    testResults$nmess[i] <- sum(testResults$nmess, na.rm=T)
    testResults[[11]] <- "messages"
    # write out to a file "test-summary.txt"
    sink(testResultsFile)
    cat(do.call("paste", lapply(testResults, format)), sep="\n")
    sink()
    # and write the same thing to "test-summary.fail" if there are any errors
    testResultsFile <- "test-summary.fail"
    if (sum(testResults$nerr, na.rm=T) > 0) {
        sink(testResultsFile)
        cat(do.call("paste", lapply(testResults, format)), sep="\n")
        sink()
    }
    return(res)
}

