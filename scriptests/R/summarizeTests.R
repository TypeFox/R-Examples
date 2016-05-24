summarizeTests <- function(debug=FALSE) {
    testResultsFile <- "test-summary.txt"
    if (file.exists(testResultsFile)) {
        all.res <- readLines(testResultsFile, -1)
        writeLines(all.res, paste(testResultsFile, ".bak", sep=""))
    } else {
        all.res <- character(0)
    }
    # all.res <- sapply(list.files(pattern=".*Rt\\.sum$"), readLines, 1)
    i <- regexpr("^.*: +[0-9]+ tests with [0-9]+ errors, [0-9]+ warnings and [0-9]+ messages", all.res) >= 1
    if (!all(i)) {
        problem.lines <- all.res[!i]
        problem.lines <- gsub(":.*$", ": ...", problem.lines, perl=TRUE)
        problem.lines <- substring(problem.lines, 1, 40)
        warning("malformed lines in '", testResultsFile, "': ", paste('"', problem.lines, '"', sep="", collapse=", "))
    }
    all.res <- all.res[i]
    all.res.con <- textConnection(all.res)
    testResults <- scan(all.res.con, quiet=T, what=list(name="", ntests=0, NULL, NULL, nerr=0, NULL, nwarn=0, NULL, NULL, nmess=0, NULL), fill=T)
    close(all.res.con)
    if (length(testResults$name)>0 && testResults$name[length(testResults$name)]=="total:")
        testResults <- lapply(testResults, "[", -length(testResults$name))
    # add a line for the totals
    i <- length(testResults$name)+1
    testResults$name[i] <- "total:"
    testResults$ntests[i] <- sum(testResults$ntests, na.rm=T)
    testResults[[3]] <- "tests"
    testResults[[4]] <- "with"
    totalErrors <- testResults$nerr[i] <- sum(testResults$nerr, na.rm=T)
    testResults[[6]] <- "errors,"
    testResults$nwarn[i] <- sum(testResults$nwarn, na.rm=T)
    testResults[[8]] <- "warnings"
    testResults[[9]] <- "and"
    testResults$nmess[i] <- sum(testResults$nmess, na.rm=T)
    testResults[[11]] <- "messages"
    i <- order(testResults$name=="total:", testResults$nerr, testResults$nwarn, testResults$nmess, testResults$name)
    testResults <- lapply(testResults, function(x, i) if (length(x)>1) x[i] else x, i)
    # write out to a file "test-summary.txt"
    # rewrite the whole output sorted with columns lined up
    sink(testResultsFile, append=FALSE)
    # cat(paste(sapply(testResults, function(x) x[length(x)]), collapse=" "), sep="\n")
    cat(do.call("paste", lapply(testResults, format)), sep="\n") # the whole thing
    sink()
    # and write the same thing to "test-summary.fail" if there are any errors
    testResultsFile <- "test-summary.fail"
    if (totalErrors > 0) {
        sink(testResultsFile)
        cat(do.call("paste", lapply(testResults, format)), sep="\n")
        sink()
        assign("test-summary.fail", envir=test.status.env, testResults)
    } else if (exists("test-summary.fail", envir=test.status.env)) {
        remove(list="test-summary.fail", envir=test.status.env)
    }
    lines <- do.call("paste", lapply(testResults, format))
    lines <- c(lines[-length(lines)], "### Overall", lines[length(lines)])
    if (totalErrors > 0) {
        firstError <- which(testResults$nerr > 0)[1]
        nFilesWithErrors <- length(lines) - firstError - 1
        lines <- c(lines[seq(1, len=firstError-1)],
                   paste("### ", nFilesWithErrors, " file", (if (nFilesWithErrors!=1) "s"),
                         " with ", totalErrors, " errors", sep=""),
                   lines[seq(firstError, length(lines))])
    } else {
        firstError <- length(testResults$nerr) # not really, but the right numbers will be output downstream
    }
    cat(paste("### Test Summary: ", firstError-1, " file", (if (firstError!=2) "s")," without errors", sep=""), lines, sep="\n")
    return(totalErrors)
}

# Environment to store status of a complete test run.
# This used to be stored in object '.test-summary.fail' in the global env,
# but that was bad practice that R CMD check --as-cran flagged.
test.status.env <- new.env()
