if (Sys.getenv('_R_CHECK_TIMINGS_')!='') {
    cat('Not running tests because appear to be running on CRAN.\n')
} else if (Sys.getenv('R_RUN_FULL_TESTS')!='1') {
    cat('Not running tests because env var R_RUN_FULL_TESTS != 1\n')
} else if (require(scriptests)) {
    runScripTests()
} else {
    cat("Not running tests because cannot find 'scriptests' package.\n")
}
