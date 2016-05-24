if (require("RUnit", quietly = TRUE)) {
    localTesting <- TRUE
    require("NMOF")
    if (localTesting)
        path <- "~/Packages/NMOF/inst/unitTests" else
    path <- system.file("unitTests", package = "NMOF")

    myTestSuite <- defineTestSuite("NMOF",
                                   dirs = path,
                                   testFileRegexp = "unitTests.+")
    stopifnot(isValidTestSuite(myTestSuite))
    testResult <- runTestSuite(myTestSuite, verbose = 0L)
    printTextProtocol(testResult, showDetails = TRUE,
                      fileName = paste(file.path(path, "report"),
                      ".txt", sep = ""))
}
