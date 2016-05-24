if(require("RUnit", quietly = TRUE))
{
    ## Coming from make test???
    if(Sys.getenv("RCMDCHECK") == "FALSE") {
        level <- Sys.getenv("LEVEL")
        pkg <- Sys.getenv("PKG")
    } else {
        ## Coming from R CMD check???
        ## level and pkg should be set if sourcing runTests.R
        ## from the command line
        if (!interactive()){
            pkg <- sub("\\.Rcheck$", '', basename(dirname(wd)))
            level <- 1
        }
    }

    library(package=pkg, character.only = TRUE)
    if(!(exists("path") && file.exists(path)))
        path <- system.file("unitTests", package = pkg)

    ## --- Testing ---

    ## Define tests
    ## level 1 are standard (default) tests run under R CMD check
    ## Add more test levels as required
    if(level == 1){
        testSuite <- defineTestSuite(name = paste(pkg, "unit testing"),
                                     dirs = path)
    } else if(level == 2){
        testSuite <- defineTestSuite(name = paste(pkg, "level 2 testing"),
                                     testFuncRegexp = "^level2test.+",
                                     dirs = path)
    } else if(level == 3){
        testSuite <- defineTestSuite(name = paste(pkg, "level 3 testing"),
                                     testFuncRegexp = "^level3test.+",
                                     dirs = path)
    } else if(level == "graphics"){
        testSuite <- defineTestSuite(name = paste(pkg, "graphics testing"),
                                     testFuncRegexp = "^graphicstest.+",
                                     dirs = path)
    }
    if(interactive()) {
        cat("Now have RUnit Test Suite 'testSuite' for package '",
            pkg, "' :\n", sep='')
        str(testSuite)
        cat('', "Consider doing",
            "\t  tests <- runTestSuite(testSuite)", "\nand later",
            "\t  printTextProtocol(tests)", '', sep = "\n")
    } else {
        ## run from shell / Rscript / R CMD Batch / ...

        if(file.access(path, 02) != 0) {
            ## cannot write to path -> use writable one
            tdir <- tempfile(paste(pkg, "unitTests", sep="_"))
            dir.create(tdir)
            pathReport <- file.path(tdir, "report")
            cat("RUnit reports are written into ", tdir, "/report.(txt|html)",
                sep = "")
        } else {
            pathReport <- file.path(path, "report")
        }

        if(level == "graphics"){
          pdf(file = paste(pathReport, "Graphics.pdf", sep = ""))
        }
        ## Run
        tests <- runTestSuite(testSuite)

        if(level == "graphics"){
          dev.off()
        }
        print(pathReport)

        ## Print Results:
        printTextProtocol(tests, showDetails = FALSE)
        printTextProtocol(tests, showDetails = FALSE,
                          fileName = paste(pathReport, "Summary.txt", sep = ""))
        printTextProtocol(tests, showDetails = TRUE,
                          fileName = paste(pathReport, ".txt", sep = ""))

        ## Print HTML Version to a File:
        ## printHTMLProtocol has problems on Mac OS X
        if (Sys.info()["sysname"] != "Darwin")
        printHTMLProtocol(tests,
                          fileName = paste(pathReport, ".html", sep = ""))

        ## stop() if there are any failures i.e. FALSE to unit test.
        ## This will cause R CMD check to return error and stop
        tmp <- getErrors(tests)
        if(tmp$nFail > 0 | tmp$nErr > 0) {
            stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
                       ", R errors: ",  tmp$nErr, ")\n\n", sep=""))
        }
    }
} else {
    cat("R package 'RUnit' cannot be loaded -- no unit tests run\n",
        "for package", pkg,"\n")
}


################################################################################
