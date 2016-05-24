
# setOldClass("RtTestSetResults")
# setOldClass("RtTestSetResultsSummary")

print.RtTestSetResultsList <- function(x, ..., transcript=FALSE) {
    if (transcript) {
        stop("can only print a transcript for one file in the list, e.g., use print(", paste(deparse(substitute(x)), collapse=" "), "[[1]], transcript=TRUE)")
    }
    print(summary(x, ...))
}

summary.RtTestSetResultsList <- function(object, ...) {
    smy <- do.call("rbind", lapply(object, function(x) {xs <- summary(x) ; n <- xs$n; data.frame(n, "tests with", errors=xs$counts[4], "errors,", warnings=xs$counts[3], "warnings and", messages=xs$counts[2], "messages")}))
    smy <- smy[order(smy$errors, smy$warnings, smy$messages), , drop=FALSE]
    smy <- rbind(smy, total=smy[1,,drop=F])
    n <- nrow(smy)
    smy[n,"n"] <- sum(smy$n[-n])
    smy[n,"errors"] <- sum(smy$errors[-n])
    smy[n,"warnings"] <- sum(smy$warnings[-n])
    smy[n,"messages"] <- sum(smy$messages[-n])
    class(smy) <- c("RtTestSetResultsList.summary", "data.frame")
    attr(smy, "dir") <- attr(object, "dir")
    attr(smy, "pattern") <- attr(object, "pattern")
    return(smy)
}

print.RtTestSetResultsList.summary <- function(x, ...) {
    cat("+++++ Test summary")
    if (!is.null(attr(x, "dir")))
        cat(" for tests in", attr(x, "dir"))
    if (!is.null(attr(x, "pattern")))
        cat("/", attr(x, "pattern"), sep="")
    cat(" +++++\n")
    cat(apply(format(cbind(paste(basename(rownames(x)), ":", sep=""), x), justify="left"), 1, paste, collapse=" "), sep="\n")
}

print.RtTestSetResults <- function(x, ..., transcript=FALSE, details=FALSE) {
    if (transcript) {
        cat("* Transcript of actual output from running commands in '", attr(x, "testname"), "':\n", sep="")
        lapply(x, function(res) {
            if (length(res$comment))
                cat(res$comment, sep="\n")
            if (length(res$transcript))
                cat(res$transcript, sep="\n")
        })
    } else if (details) {
        lapply(x, function(res) {
            if (res$status=="ok")
                cat(".")
            else
                cat("\n", paste(res$msg, "\n", sep=""), sep="")
        })
        cat("\n")
    } else {
        print(summary(x, ...))
    }
    invisible(x)
}

summary.RtTestSetResults <- function(object, ...) {
    y <- table(factor(sapply(object, "[[", "status"), levels=c("ok", "info", "warning", "error")))
    n <- object[[length(object)]]$i
    structure(list(n=n, counts=y, testname=basename(attr(object, "testname"))), class="RtTestSetResults.summary")
}

print.RtTestSetResults.summary <- function(x, ...) {
    cat(x$testname, ": ", x$n, " tests with ",
        x$counts["error"], " errors, ",
        x$counts["warning"], " warnings and ",
        x$counts["info"], " messages\n", sep="")
    invisible(x)
}

# is setupTests() needed for anything?
# setupTests <- function(debug=FALSE, create.Rout.save=FALSE, addSelfCheck=FALSE) initializeTests(debug, create.Rout.save, addSelfCheck)


if (F) {
    # was using this in an earlier incarnation
testWrapper <- function(pkg.dir, dir=file.path(paste(pkg.dir, ".Rcheck", sep=""), "tests")) {
    cwd <- getwd()
    on.exit(setwd(cwd), add=TRUE)
    setwd(dir)
    existing.files <- list.files()
    status <- .runPackageTests()
    new.files <- setdiff(list.files(), existing.files)
    if (length(new.files)) {
        cat("* Removing ", length(new.files), " new files: ", paste(new.files, collapse=", "), "\n", sep="")
        file.remove(new.files)
    }
    return(status)
}
}

