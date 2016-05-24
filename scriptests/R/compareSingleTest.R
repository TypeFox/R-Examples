compareSingleTest <- function(input, control, target, comment, garbage, actual, testnum, testname, verbose) {
    ## Return a list containing the following elements:
    ##   status: one of 'ok', 'error', 'warn', 'info'
    ##   msg: character vector describing differences (for human readability)
    ## A zero length return value indicates a match.
    ## Strip trailing & leading white space, convert all whitespace to single, remove blank lines.
    ## If the control text contains a line starting with "#@ ignore output" (case-insensitive),
    ## then the match is successful, regardless of the actual output.
    ## All lines in the control text that start with "#@ gsub(pattern, replacement, var)",
    ## where 'var' is 'target' or 'actual' (no quotes) are applied as substitutions, in order,
    ## to the target or actual text, before attempting a match.
    mismatch.status <- "error"
    status.msg <- NULL
    ignore.whitespace <- TRUE
    ignore.linebreaks <- FALSE
    diff.msg <- NULL # message to output if there are differences
    did.gsub <- c(FALSE, FALSE)
    if (!is.null(control)) {
        if (any(regexpr('^\\#@ *ignore[- ]output', control, ignore.case=TRUE)>0))
            target <- actual <- character(0)
        for (line in control) {
            if (regexpr('^\\#@ *gsub\\(', casefold(line, upper=FALSE))>0) {
                expr <- try(parse(text=gsub("^\\#@ *", "", line, perl=TRUE)), silent=TRUE)
                v <- all.vars(expr)
                if (is(expr, "try-error")) {
                    warning("parse error in test-control line '", line, "': ",
                            gsub("Error in parse\\(file, n, text, prompt\\) :[ \n\t]+", "", expr, perl=TRUE))
                } else if (!all(is.element(v, c("target", "actual", "both")))) {
                    warning("ignoring test-control line '", line, "' uses variables other than 'target', 'actual' and 'both'")
                } else if (length(v)==0) {
                    warning("ignoring test-control line '", line, "' that refers to neither 'target', 'actual' nor 'both'")
                } else if (length(v)>1) {
                    warning("ignoring test-control line '", line, "' that refers to more than one of 'target', 'actual' and 'both'")
                } else if (v=="target") {
                    res <- try(eval(do.call("substitute", list(expr[[1]], list(target=target)))), silent=TRUE)
                    if (is(res, "try-error"))
                        warning("error in running gsub control line '", line, "': ", as.character(res))
                    else
                        target <- res
                    did.gsub[1] <- TRUE
                } else if (v=="actual") {
                    res <- try(eval(do.call("substitute", list(expr[[1]], list(actual=actual)))), silent=TRUE)
                    if (is(res, "try-error"))
                        warning("error in running gsub control line '", line, "': ", as.character(res))
                    else
                        actual <- res
                    did.gsub[2] <- TRUE
                } else if (v=="both") {
                    res <- try(eval(do.call("substitute", list(expr[[1]], list(both=target)))), silent=TRUE)
                    if (is(res, "try-error"))
                        warning("error in running gsub control line '", line, "': ", as.character(res))
                    else
                        target <- res
                    res <- try(eval(do.call("substitute", list(expr[[1]], list(both=actual)))), silent=TRUE)
                    if (is(res, "try-error"))
                        warning("error in running gsub control line '", line, "': ", as.character(res))
                    else
                        actual <- res
                    did.gsub[1:2] <- TRUE
                }
            } else if (regexpr('^\\#@ *warn[- ]only', line, ignore.case=TRUE)>0) {
                mismatch.status <- "warning"
                msg <- gsub('^\\#@ *warn[- ]only:? ?', '', line, perl=TRUE)
                if (msg != line && nchar(line)>0)
                    diff.msg <- msg
                status.msg <- c(status.msg, line)
            } else if (regexpr('^\\#@ *info[- ]only', line, ignore.case=TRUE)>0) {
                mismatch.status <- "info"
                msg <- gsub('^\\#@ *info[- ]only:? ?', '', line, perl=TRUE)
                if (msg != line && nchar(line)>0)
                    diff.msg <- msg
                status.msg <- c(status.msg, line)
            } else if (regexpr('^\\#@ *keep[- ]whitespace', line, ignore.case=TRUE)>0) {
                ignore.whitespace <- FALSE
            } else if (regexpr('^\\#@ *ignore[- ]linebreaks', line, ignore.case=TRUE)>0) {
                ignore.linebreaks <- TRUE
            } else if (regexpr('^\\#@ *diff[- ]msg:', line, ignore.case=TRUE)>0) {
                diff.msg <- c(diff.msg, gsub('^\\#@ *diff-msg:? ?', '', line, perl=TRUE))
            } else if (regexpr('^\\#@ *ignore[- ]output', line, ignore.case=TRUE)>0) {
            } else if (regexpr('^\\#@#', line)<1) {
                warning("cannot understand test-control line (ignoring): '", line, "'")
            }
        }
    }
    if (ignore.linebreaks) {
        ## Remove newlines and join lines if necessary.
        target <- gsub("\n", "", paste(target, collapse = ""), perl=TRUE)
        actual <- gsub("\n", "", paste(actual, collapse = ""), perl=TRUE)
    }
    # cannonicalize whitespace
    if (ignore.whitespace) {
        target <- gsub("^[[:space:]]*", "", target, perl=TRUE)
        target <- gsub("[[:space:]]*$", "", target, perl=TRUE)
        target <- gsub("[[:space:]]{2,}", " ", target, perl=TRUE)
        target <- target[target!=""]
        actual <- gsub("^[[:space:]]*", "", actual, perl=TRUE)
        actual <- gsub("[[:space:]]*$", "", actual, perl=TRUE)
        actual <- gsub("[[:space:]]{2,}", " ", actual, perl=TRUE)
        actual <- actual[actual!=""]
    }
    target.len <- length(target)
    actual.len <- length(actual)
    if (length(target) < length(actual))
        target <- c(target, rep("", length(actual) - length(target)))
    if (length(target) > length(actual))
        actual <- c(actual, rep("", length(target) - length(actual)))
    i <- target == actual
    if (length(i) && !all(i)) {
        msg <- c(paste("* ", casefold(substring(mismatch.status, 1, 1), upper=TRUE), substring(mismatch.status, 2),
                       " mismatch on output for test number ", testnum, " in ", testname, ":", sep=""),
                 status.msg, input)
        if (!is.null(diff.msg))
            msg <- c(msg, diff.msg)
        if (target.len==0) {
            msg <- c(msg, paste("* Empty target output, actual output is:"), paste("", actual, sep=""))
        } else if (actual.len==0) {
            msg <- c(msg, paste("* No actual output, but target output", if (length(target)>3) "began with:" else "is:"),
                     target[seq(len=min(3, length(target)))])
        } else if (actual.len+target.len <= 10) {
            msg <- c(msg, (if (did.gsub[1]) "* Target output (after gsub):" else "* Target output:"),
                     paste("", target[seq(target.len)], sep=""),
                     (if (did.gsub[2]) "* Actual output (after gsub):" else "* Actual output:"),
                     paste("", actual[seq(actual.len)], sep=""))
        } else {
            msg <- c(msg, paste("* Output differs first at line ", which(!i)[1], ":", sep=""),
                     paste((if (did.gsub[1]) "  target (after gsub):" else "  target:"), target[which(!i)[1]]),
                     paste((if (did.gsub[2]) "  actual (after gsub):" else "  actual:"), actual[which(!i)[1]]))
        }
        status <- mismatch.status
    } else {
        msg <- character(0)
        status <- "ok"
    }
    if (length(garbage)) {
        if (is.null(input)) {
            msg <- c(msg, paste("* Warning:", length(garbage), "uninterpretable line(s) around comment",
                                testnum, "in", testname, ":"), comment,
                     "* Uninterpretable lines are:", garbage)
        } else {
            if (status=="ok")
                msg <- c(msg, paste("* Warning:", length(garbage), "uninterpretable line(s) around test number",
                                    testnum, "in", testname, ":"), input,
                         "* Uninterpretable lines are:", garbage)
            else
                msg <- c(msg, paste("* Warning:", length(garbage), "uninterpretable line(s) around the above test case:"), garbage)
        }
        if (status=="ok" || status=="info")
            status <- "warning"
    }
    if (verbose)
        if (status=="ok")
            cat(".")
        else
            cat("", msg, sep="\n")
    res <- list(status=status, msg=msg, i=testnum)
    return(res)
}

