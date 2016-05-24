context("FitAnalyzer tests")

test_that("FitAnalyzer test cases", {
    
    masterPath <-
        system.file("extdata", "daily-time-series", package = "fitcoach")

    ### Tests for daily analysis
  
    # Test 1 - Initializing and dataframe creation
    ana <- FitAnalyzer$new("calories")
    ts <-
        ana$getAnalysisFrame(folder = masterPath, analysis.type = "daily")
    expect_equal(nrow(ts), 191)
    
    # Test 2 - Find important variables
    vars <- ana$findImportantVariables(tsDataFrame = ts, seed = 12345)
    expect_equal(vars$name[1], "minutesLightlyActive")
    
    # Test 3 - Find important variables, without arguments
    vars <- ana$findImportantVariables()
    expect_equal(names(vars[1]), "Overall")
    
    # Test 4 - Goal prediction
    rows.test <- ts[c(3,7,10), ]
    rows.test <- createDependentVariableFrame(master = rows.test, goal = "calories")
    res <- ana$predictGoal(rows.test)
    expect_less_than(res[1], 2820)
    expect_gte(res[1], 2819)
    
    ### Tests for intra-day analysis
    
    # Test 5 - Initializing and dataframe creation
    masterPath <-
        system.file("extdata", "intra-daily-timeseries", package = "fitcoach")
    ana <- FitAnalyzer$new("calories")
    intra <-
        ana$getAnalysisFrame(folder = masterPath, analysis.type = "intra.day")
    expect_equal(nrow(intra), 2016)

    # Test 6 - Find important variables
    vars <- ana$findImportantVariables(intra)
    vars <- sort(vars, decreasing = TRUE)
    expect_equal(names(vars[1]), "steps")
    
    # Test 7 - Predict goal for the day
    rows.test <- intra[c(3), ]
    res <- ana$predictGoal(rows.test)
    expect_less_than(res , 2518)

})
