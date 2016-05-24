initializeTests <- function(debug=FALSE, create.Rout.save=FALSE, addSelfCheck=FALSE, pattern=NULL, subst=NULL, R.suf="R") {
    R.suf.regexp <- paste("\\.", R.suf, "$", sep="")
    # Create .R and .Rout.save files for each .Rt file
    wd <- getwd()
    if (debug)
        cat("   initializeTests: wd=", wd, "\n")
    # Expect that we are running in .../<pkg.dir>.Rcheck/tests
    test.dir <- basename(dirname(wd)) # should get something like mypackage.Rcheck
    pkg.name <- NULL
    if (debug)
        cat("   initializeTests: debug=", debug, "\n")
    if (regexpr("\\.Rcheck$", test.dir) > 0) {
        # Can't rely on "*" being the package name in "*.Rcheck", e.g.,
        # in r-forge directory structure, packages are structured like this:
        # mypackage/pkg/{DESCRIPTION,R,man} etc.
        # Depending on how tests are run, tests go into mypackage/pkg.Rcheck
        # or into mypackage/mypackage.Rcheck.
        # So get the package name from the log messages in <pkg>.Rcheck/00install.out
        if (debug)
            cat("   Looking for ", file.path(dirname(wd), "00install.out"), "\n", sep="")
        if (file.exists(file.path(dirname(wd), "00install.out"))) {
            pkg.name <- gsub(")", "", gsub("* DONE (", "", fixed=TRUE, grep("* DONE ", readLines(file.path(dirname(wd), "00install.out")), fixed=TRUE, value=TRUE)))
            if (length(pkg.name)>1)
                pkg.name < pkg.name[nchar(pkg.name)>0]
            if (length(pkg.name)>1) {
                warning("found several lines matching '* DONE (...)' in ", file.path(dirname(wd), "00install.out"))
                pkg.name <- pkg.name[1]
            }
            if (debug)
                cat("   Read pkg.name= '", pkg.name, "'\n", sep="")
            if (length(pkg.name)<1)
                cat("   Failed to work out pkg.name from 00install.line: ", grep("DONE", readLines(file.path(dirname(wd), "00install.out")), fixed=TRUE, value=TRUE), "\n")
        }
        if (length(pkg.name) && nchar(pkg.name)==0)
            pkg.name <- NULL
    }

    # process all .R files - if a .Rout.save file exists, generate a .Rt.save file
    for (cmdIn in list.files(pattern=R.suf.regexp, ignore.case=TRUE)) {
        if (!is.null(pattern) && !length(grep(pattern, cmdIn))) {
            if (debug)
                cat(" Skipping file ", cmdIn, "\n")
            next
        }
        rOutSave <- gsub(R.suf.regexp, ".Rout.save", cmdIn, perl=TRUE)
        if (file.exists(rOutSave)) {
            rtSave <- gsub(R.suf.regexp, ".Rt.save", cmdIn, perl=TRUE)
            if (debug)
                cat(" Pre-processing tests in ", cmdIn, "/", rOutSave, " to generate ", rtSave, "\n", sep="")
            tests <- parseTranscriptFile(rOutSave)
            env <- new.env()
            save(list="tests", file=rtSave, envir=env)
        } else {
            if (debug)
                cat(" No saved output file for tests in ", cmdIn, "; will just check that tests run without stopping with an error\n")
        }
    }

    # process all .Rt files (generate a .R and .Rt.save and optionally .Rout.save)
    for (rtIn in list.files(pattern=".*\\.Rt$", ignore.case=TRUE)) {
        if (!is.null(pattern) && !length(grep(pattern, rtIn))) {
            if (debug)
                cat(" Skipping file ", rtIn, "\n")
            next
        }
        rOutSave <- gsub("\\.Rt$", ".Rout.save", rtIn, perl=TRUE)
        rtSave <- gsub("\\.Rt$", ".Rt.save", rtIn, perl=TRUE)
        cmdOut <- gsub("\\.Rt$", paste(".", R.suf, sep=""), rtIn, perl=TRUE)
        if (debug)
            cat(" Pre-processing tests in ", rtIn, " to generate ", cmdOut, if (create.Rout.save) paste(",", rOutSave), ", and ", rtSave, "\n", sep="")
        tests <- parseTranscriptFile(rtIn, subst=subst)
        env <- new.env()
        test.obj.name <- paste(gsub("\\.Rt$", "", rtIn), ".tests", sep="", perl=TRUE)
        assign(test.obj.name, tests, envir=env)
        save(list="tests", file=rtSave, envir=env)
        ## add a comment to the beginning for what is printed in the file
        tests <- c(list(list(comment=paste("> # generated automatically in", getwd(), "on", Sys.time())),
                        if (length(pkg.name)) list(input=paste("> library('", pkg.name, "', char=TRUE)", sep=""))
                        else list(input=paste("> # could not work out package name from getwd: '", dirname(wd), "/00install.out'", sep="")),
                        list(input="> searchpaths() # seeing where these came from can be useful for debugging"),
                        list(comment="> # End of scriptests preamble")),
                   tests)

        ## write out the commands and the desired output
        cmdOutCon <- file(cmdOut, "w")
        if (create.Rout.save)
            rOutSaveCon <- file(rOutSave, "w")
        # do in a local() block so we can use on.exit()
        local({
            on.exit(close(cmdOutCon), add=TRUE)
            if (create.Rout.save)
                on.exit(close(rOutSaveCon), add=TRUE)
            lapply(tests, function(test) {
                if (!is.null(test$comment)) {
                    commands <- test$comment
                    commandsIn <- gsub("^[>+] ?", "", grep("^[>+]", commands, value=TRUE, perl=TRUE), perl=TRUE)
                    output <- NULL
                } else if (!is.null(test$input)) {
                    commands <- test$input
                    output <- test$output
                    commandsIn <- gsub("^[>+] ?", "", grep("^[>+]", commands, value=TRUE, perl=TRUE), perl=TRUE)
                } else if (!is.null(test$garbage)) {
                    commands <- test$garbage
                    commandsIn <- character(0)
                    output <- NULL
                    return(NULL)
                } else {
                    ## stop("invalid test structure: no 'comment', 'input', or 'garbage' component")
                    return(NULL)
                }
                cat(file=cmdOutCon, paste(commandsIn, "\n", sep=""), sep="")
                if (create.Rout.save)
                    cat(file=rOutSaveCon, paste(commands, "\n", sep=""), sep="")
                if (create.Rout.save && length(output))
                    cat(file=rOutSaveCon, paste(output, "\n", sep=""), sep="")

            })
            if (addSelfCheck) {
                cat(file=cmdOutCon, "# End of scriptests output\n")
                cat(file=cmdOutCon, "flush(stdout())\n")
                cat(file=cmdOutCon, "require('scriptests', character.only=TRUE)\n")
                cat(file=cmdOutCon, "scriptests:::checkTestOutput('", rtIn, "', '", rtSave, "')\n", sep="")
            }
            close(cmdOutCon)
            if (create.Rout.save)
                close(rOutSaveCon)
            on.exit()
        })
    }
    return(0)
}
