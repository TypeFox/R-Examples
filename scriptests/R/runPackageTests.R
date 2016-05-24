## Modification of function .runPackageTests() from R-2.9.0/src/library/tools/R/testing.R
## used by R CMD check
## This is a private function in scriptests, but it was originally written so that
## it could be a drop in replacement for .runPackageTests() in .../src/library/tools/R/testing.R
##
## run.from is the name of the file from which the tests were called -- don't
## want to run it again because that way lies infinite recursion!

.runPackageTests <- function(use_gct = FALSE, run.preexisting.R.files=TRUE,
                             initializeFun=NULL, finalizeFun=NULL, diffFun=NULL,
                             debug=FALSE, stopOnError=TRUE, pattern=NULL, subst=NULL,
                             run.from=NULL)
{
    runone <- function(f, diffFun=NULL, stopOnError=TRUE, debug=TRUE, R.suf="Rs")
    {
        cat("** Running ", sQuote(f), " in ", getwd(), "\n", sep="")
        R.suf.regexp <- paste("\\.", R.suf, "$", sep="")
        outfile <- gsub(R.suf.regexp, ".Rout", f)
        cmd <- paste(shQuote(file.path(R.home(), "bin", "R")),
                     "CMD BATCH --vanilla --no-timing",
                     shQuote(f), shQuote(outfile))
        if (.Platform$OS.type == "windows") {
            Sys.setenv(LANGUAGE="C")
            Sys.setenv(R_TESTS="startup.Rs")
        } else {
            cmd <- paste("LANGUAGE=C", "R_TESTS=startup.Rs", cmd)
        }
        if (debug)
            cat("   Running cmd: ", cmd, "\n", sep="")
        res <- system(cmd)
        if (res) {
            cat("   Failure running ", sQuote(f), ": returned code ", res, "\n")
            if (stopOnError) {
                file.rename(outfile, paste(outfile, "fail", sep="."))
                return(1L)
            }
        }
        # Look for a .Rout.save file -- this can be supplied by the user
        savefile <- paste(outfile, "save", sep = "." )
        # Look for a .Rt.save file -- this was auto-generated and we want to remove it
        # after we've used it.
        rtSave <- gsub(R.suf.regexp, ".Rt.save", f, perl=TRUE)
        if (file.exists(rtSave)) {
            on.exit(unlink(rtSave))
        }
        if (is.null(diffFun)) {
            if (file.exists(savefile)) {
                cat("   Comparing ", sQuote(outfile), " to ", sQuote(savefile), " ...")
                res <- Rdiff(outfile, savefile, TRUE)
                if (!res) cat(" OK\n")
            }
        } else {
            args <- list(commandfile=f, outfile=outfile)
            if (file.exists(savefile))
                args <- c(args, list(savefile=savefile))
            for (i in seq(to=2, by=-1, len=max(0, length(diffFun)-1))) {
                diffFun[i+length(args)] <- diffFun[i]
                names(diffFun)[i+length(args)] <- names(diffFun[i])
            }
            diffFun[seq(2,len=length(args))] <- args
            names(diffFun)[seq(2,len=length(args))] <- names(args)
            if (debug)
                cat("   Calling ", format(diffFun), "\n", sep="")
            res <- try(eval(diffFun))
            if (inherits(res, "try-error")) {
                cat("   Error evaluting ", format(diffFun), ": ", as.character(res), "\n", sep="")
                return(1L)
            } else if (!is.numeric(res)) {
                cat("   Error: ", format(diffFun), " returned non-numeric result:", as.character(res), "\n", sep="")
                return(1L)
            }
        }
        return(0L)
    }
    if (run.preexisting.R.files) {
        if (Sys.getenv("TESTS_HAVE_RUN_R_FILES")=="") {
            Sys.setenv("TESTS_HAVE_RUN_R_FILES", "TRUE")
        } else {
            warning("refusing to run pre-existing .R and .Rin files: env var TESTS_HAVE_RUN_R_FILES is non-empty, so running these files might cause an infinite regression")
            run.preexisting.R.files <- FALSE
        }
    }
    # Can change the default R command file suffix here -- a reason to change
    # it would be that if the creation of .R files in the tests directory is
    # confusing other aspects of the testing.  Setting R.suf to a different
    # value will make scriptests use that value as the suffix for R command files.
    # R.suf should NOT begin with a dot.
    R.suf <- "R"
    preexisting.Rin.files <- dir(".", pattern="\\.Rin$", ignore.case=TRUE)
    preexisting.R.files <- dir(".", pattern="\\.R$", ignore.case=TRUE)
    preexisting.R.files <- setdiff(preexisting.R.files, run.from)

    file.copy(file.path(R.home("share"), "R", "tests-startup.R"), "startup.Rs")
    if (use_gct) cat("gctorture(TRUE)" , file = "startup.Rs", append = TRUE)
    # Need to make consistent behavior between interactive and non-interactive
    # for options(showErrorCalls), and can't get evalCapture() to capture the
    # error calls, so turn them off for non-interactive.
    cat("options(showErrorCalls = FALSE)\n", file = "startup.Rs", append = TRUE)
    nfail <- 0L ## allow for later running all tests even if some fail.

    needPkg <- character(0)
    if (file.exists("CONFIG")) {
        # need to make sure that the "tools" package is attached
        config <- read.dcf("CONFIG")[1L,]
        names(config) <- casefold(names(config), upper=FALSE)

        if (is.element("debug", names(config))) {
            debug <- try(eval(parse(text=config["debug"])[[1]])!=0)
            if (inherits(debug, "try-error")) {
                warning("could not interpret 'Debug' entry in CONFIG (\"", config["debug"], "\") as a logical value")
                debug <- TRUE
            } else if (length(debug)!=1 || is.na(debug)) {
                warning("could not interpret 'Debug' entry in CONFIG (\"", config["debug"], "\") did not evaluate to a single TRUE/FALSE value")
                debug <- TRUE
            }
            if (debug)
                cat("   Setting debug=TRUE\n")
        }

        # a line like "Config: x:::y" in the CONFIG file will call x:::y to get config info
        if (is.element("config", names(config))) {
            configFun <- try(parse(text=config["config"])[[1]])
            configPkg <- character(0)
            if (!is.call(configFun)) {
                warning("could not parse 'Config' entry in CONFIG (\"", config["config"], "\") as a function call")
                return(1)
            } else if (!is.name(configFun[[1]]) && as.character(configFun[[1]][[1]]) == ":::") {
                configPkg <- as.character(configFun[[1]][[2]])
            } else {
                warning("'Config' entry in CONFIG (\"", config["config"], "\") must be function call like 'package:::fun()'")
                return(1)
            }
            for (pkg.name in setdiff(configPkg, .packages())) {
                if (debug)
                    cat("Attemping to load package '", pkg.name, "'\n")
                if (!library(pkg.name, character.only=TRUE, logical.return=TRUE)) {
                    warning("could not load package '", pkg.name, "' needed for testing initialize/diff/finalize calls")
                    return(1)
                }
            }
            cat("   Using Config =", format(configFun), "\n")
            if (debug)
                cat("   Calling config function ", format(configFun), "\n", sep="")
            res <- try(eval(configFun), silent=TRUE)
            if (inherits(res, "try-error")) {
                warning("failed to run config function call ", format(configFun), ": ", paste(as.character(res), collapse=" "))
                return(1)
            }
            if (!is.character(res) || is.null(names(res))) {
                warning("Config function ", format(configFun), " did not return a named character vector (returned ", paste(as.character(res), collapse=" "), ")")
                return(1)
            }
            names(res) <- casefold(names(res), upper=FALSE)
            config <- c(config, res)
        }

        depends <- character(0)
        # should we work out what the name of the package being tested, and add that too?
        if (is.element('Depends', names(config))) {
            # old code: depends <- tools:::.get_requires_from_package_db(config, "Depends")
            depends <- unlist(strsplit(config['Depends'], ","))
            depends <- setdiff(sub("^[[:space:]]*([[:alnum:].]+).*$", "\\1", depends), 'R')
        }

        if (is.element("stoponerror", names(config))) {
            stopOnError <- try(eval(parse(text=config["stoponerror"])[[1]])!=0)
            if (inherits(stopOnError, "try-error")) {
                warning("could not interpret 'StopOnError' entry in CONFIG (\"", config["stoponerror"], "\") as a logical value")
                stopOnError <- FALSE
            } else if (length(stopOnError)!=1 || is.na(stopOnError)) {
                warning("could not interpret 'StopOnError' entry in CONFIG (\"", config["stoponerror"], "\") did not evaluate to a single TRUE/FALSE value")
                stopOnError <- FALSE
            } else {
                cat("   Setting stopOnError=", stopOnError, "\n")
            }
        }

        if (length(depends))
            cat(paste("library('", unique(depends), "', character.only=TRUE)\n", sep=""), file = "startup.Rs", append = TRUE)

        if (is.element("initialize", names(config))) {
            initializeFun <- try(parse(text=config["initialize"])[[1]])
            if (!is.call(initializeFun)) {
                warning("could not parse 'Initialize' entry in CONFIG (\"", config["initialize"], "\") as a function call")
                return(1)
            } else {
                if (!is.name(initializeFun[[1]]) && as.character(initializeFun[[1]][[1]]) == ":::")
                    needPkg <- unique(c(needPkg, as.character(initializeFun[[1]][[2]])))
                cat("   Using Initialize =", format(initializeFun), "\n")
            }
        }
        if (is.element("finalize", names(config))) {
            finalizeFun <- try(parse(text=config["finalize"])[[1]])
            if (!is.call(finalizeFun)) {
                warning("could not parse 'Finalize' entry in CONFIG (\"", config["finalize"], "\") as a function call")
                return(1)
            } else {
                if (!is.name(finalizeFun[[1]]) && as.character(finalizeFun[[1]][[1]]) == ":::")
                    needPkg <- unique(c(needPkg, as.character(finalizeFun[[1]][[2]])))
                cat("   Using Finalize =", format(finalizeFun), "\n")
            }

        }
        if (is.element("diff", names(config))) {
            diffFun <- try(parse(text=config["diff"])[[1]])
            if (!is.call(diffFun)) {
                warning("could not parse 'Diff' entry in CONFIG (\"", config["diff"], "\") as a function call")
                return(1)
            } else {
                if (!is.name(diffFun[[1]]) && as.character(diffFun[[1]][[1]]) == ":::")
                    needPkg <- unique(c(needPkg, as.character(diffFun[[1]][[2]])))
                cat("   Using Diff =", format(diffFun), "\n")
            }
        }
        if (is.element("rsuffix", names(config))) {
            R.suf <- config["rsuffix"]
            if (inherits(R.suf, "try-error")) {
                warning("could not interpret 'Rsuffix' entry in CONFIG (\"", config["rsuffix"], "\") as a logical value")
                R.suf <- "Rs"
            } else if (length(R.suf)!=1 || is.na(R.suf)) {
                warning("could not interpret 'Rsuffix' entry in CONFIG (\"", config["rsuffix"], "\") did not evaluate to a single character value")
                R.suf <- "Rs"
            }
            if (debug)
                cat("   Setting Rsuffix=", R.suf, "\n")
        }
    }

    if (length(needPkg)) {
        for (pkg.name in setdiff(needPkg, .packages())) {
            if (debug)
                cat("Attemping to load package '", pkg.name, "'", "\n")
            if (!library(pkg.name, character.only=TRUE, logical.return=TRUE)) {
                warning("could not load package '", pkg.name, "' needed for testing initialize/diff/finalize calls")
                return(1)
            }
        }
    }

    cat("### Show output from here ###\n")
    if (!is.null(initializeFun)) {
        # Add a pattern=<pattern> actual argument to the call to initializeFun
        # if it does have a formal argument named 'pattern'
        init.args <- names(formals(eval(initializeFun[[1]])))
        if (!is.null(pattern) && is.element("pattern", init.args)) {
            initializeFun[length(initializeFun)+1] <- pattern
            names(initializeFun)[length(initializeFun)] <- "pattern"
        }
        # Add a subst=<subst> actual argument to the call to initializeFun
        # if it does have a formal argument named 'subst'
        if (!is.null(subst) && is.element("subst", init.args)) {
            initializeFun[length(initializeFun)+1] <- subst
            names(initializeFun)[length(initializeFun)] <- "subst"
        }
        # Add a debug=<debug> actual argument to the call to initializeFun
        # if it does have a formal argument named 'debug'
        if (!is.null(debug) && !identical(debug, FALSE) && is.element("debug", init.args) && !is.element("debug", names(initializeFun))) {
            initializeFun[length(initializeFun)+1] <- debug
            names(initializeFun)[length(initializeFun)] <- "debug"
        }
        # Add a R.suf=<R.suf> actual argument to the call to initializeFun
        # if it does have a formal argument named 'R.suf'
        if (!is.null(R.suf) && is.element("R.suf", init.args)) {
            initializeFun[length(initializeFun)+1] <- R.suf
            names(initializeFun)[length(initializeFun)] <- "R.suf"
        }
        if (debug)
            cat("   Calling initialization function ", format(initializeFun), "\n", sep="")
        res <- try(eval(initializeFun), silent=TRUE)
        if (inherits(res, "try-error")) {
            warning("failed to run initialize function call ", format(initializeFun), ": ", paste(as.character(res), collapse=" "))
            return(1)
        }
        if (!is.atomic(res) || length(res)!=1 || is.na(res) || res != 0) {
            warning("initialize function call ", format(initializeFun), " ran successfully but returned non-zero value ",
                    if (length(res)!=1) paste("(length=", length(res), ")", sep=""), ": ",
                    paste(format(res), collapse="; "))
            return(1)
        }
    }
    if (!is.null(diffFun)) {
        # Add a R.suf=<R.suf> actual argument to the call to diffFun
        # if it does have a formal argument named 'R.suf'
        if (!is.null(R.suf) && is.element("R.suf", init.args)) {
            diffFun[length(diffFun)+1] <- R.suf
            names(diffFun)[length(diffFun)] <- "R.suf"
        }
        # Add a debug=<debug> actual argument to the call to diffFun
        # if it does have a formal argument named 'debug'
        if (!is.null(debug) && !identical(debug, FALSE) && is.element("debug", init.args) && !is.element("debug", names(diffFun))) {
            diffFun[length(diffFun)+1] <- debug
            names(diffFun)[length(diffFun)] <- "debug"
        }
    }

    R.suf.regexp <- paste("\\.", R.suf, "$", sep="")

    Rinfiles <- dir(".", pattern="\\.Rin$", ignore.case=TRUE)
    if (!run.preexisting.R.files) {
        ## message("ignoring run.from=", run.from)
        if (debug && length(preexisting.Rin.files))
            cat("   Not running these pre-existing .Rin files: ", paste(preexisting.Rin.files, collapse=" "), "\n")
        Rinfiles <- setdiff(Rinfiles, preexisting.Rin.files)
    }
    for(f in Rinfiles) {
        Rfile <- sub("\\.Rin$", paste(".", R.suf, sep=""), f)
        cat("   Creating ", sQuote(Rfile), " from ", f, "\n")
        cmd <- paste(shQuote(file.path(R.home(), "bin", "R")),
                     "CMD BATCH --no-timing --vanilla --slave", f)
        if (system(cmd))
            warning("creation of ", sQuote(Rfile), " failed")
        else if (file.exists(Rfile)) nfail <- nfail + runone(Rfile, stopOnError=stopOnError, debug=debug, R.suf=R.suf)
        if (nfail > 0) return(nfail)
    }

    Rfiles <- dir(".", pattern=R.suf.regexp, ignore.case=TRUE)
    Rfiles <- setdiff(Rfiles, c(run.from, "startup.Rs"))
    if (!run.preexisting.R.files) {
        Rfiles <- setdiff(Rfiles, preexisting.R.files)
        if (debug && length(preexisting.R.files))
            cat("   Not running these pre-existing .R files: ", paste(preexisting.R.files, collapse=" "), "\n")
    }
    for(f in Rfiles) {
        nfail <- nfail + runone(f, diffFun, stopOnError=stopOnError, debug=debug, R.suf=R.suf)
        if (nfail > 0) return(nfail)
    }

    if (!is.null(finalizeFun)) {
        # Add a debug=<debug> actual argument to the call to finalizeFun
        # if it does have a formal argument named 'debug'
        if (!is.null(debug) && !identical(debug, FALSE) && is.element("debug", init.args) && !is.element("debug", names(finalizeFun))) {
            finalizeFun[length(finalizeFun)+1] <- debug
            names(finalizeFun)[length(finalizeFun)] <- "debug"
        }
        if (debug)
            cat("   Calling finalization function ", format(finalizeFun), "\n", sep="")
        res <- try(eval(finalizeFun), silent=TRUE)
        if (inherits(res, "try-error")) {
            warning("failed to run finalize function call ", format(finalizeFun), ": ", paste(as.character(res), collapse=" "))
            nfail <- nfail + 1
        } else if (!is.atomic(res) || length(res)!=1 || is.na(res) || res != 0) {
            nfail <- nfail + 1
            if (!file.exists("test-summary.fail")) {
                # don't need to give multiple indications of problems
                warning("finalize function call ", format(finalizeFun), " ran successfully but returned non-zero value ",
                        if (length(res)!=1) paste("(length=", length(res), ")", sep=""), ": ",
                        paste(format(res), collapse="; "))
             }
        }
    }

    return(nfail)
}

