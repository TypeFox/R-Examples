# unit testing
if (require(RUnit, quietly = TRUE)) {

    # load package
    library(spcosa)

    # define the test suite
    testSuite <- defineTestSuite(
        name = "spcosa",
        dirs = file.path(getwd(), "unitTests"),
        testFileRegexp = "^runit_",
        testFuncRegexp = "^test_"
    )
    
    # run test suite
    testResult <- runTestSuite(testSuite)
    
    # print HTML version to a file
    printHTMLProtocol(testResult, fileName = "unitTesting.html")
}
