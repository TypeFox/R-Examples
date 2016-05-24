testsuite.bootES <- defineTestSuite("bootES",
                                    dirs=system.file("tests", package="bootES"),
                                    testFileRegexp="^runit",
                                    testFuncRegexp="^test")
