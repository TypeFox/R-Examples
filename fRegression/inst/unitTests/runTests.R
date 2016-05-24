
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################


pkg <- "fRegression"

if(require("RUnit", quietly = TRUE))
{

    library(package=pkg, character.only = TRUE)
    if(!(exists("path") && file.exists(path)))
        path <- system.file("unitTests", package = pkg)

    ## --- Testing ---

    ## Define tests
    testSuite <- defineTestSuite(name = paste(pkg, "unit testing"),
                                 dirs = path)

    if(interactive()) {
        cat("Now have RUnit Test Suite 'testSuite' for package '",
            pkg, "' :\n", sep='')
        str(testSuite)
        cat('', "Consider doing",
            "\t  tests <- runTestSuite(testSuite)", "\nand later",
            "\t  printTextProtocol(tests)", '', sep = "\n")
    } else {
        ## run from shell / Rscript / R CMD Batch / ...
        ## Run
        tests <- runTestSuite(testSuite)

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
