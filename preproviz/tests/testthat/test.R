library(preproviz)
context("Object initialization")

testdata <- data.frame(matrix(rnorm(3*20, 1, .5), ncol=3),class=sample(letters[1:2], 20, replace=TRUE))
testdataobject <- initializedataobject(testdata)
testmissingvalueshare <- initializesubclassobject("MissingValueShare", testdataobject)
testobjectsfromparameters <- getinitializedsubclassobjects(testdataobject, defaultParameters)
testmissingvalueshareisvalid <- testmissingvalueshare@isvalid
#testfeature <- constructfeature("testfeature", "apply(data, 1, function(x) sum(is.na))", impute=TRUE)
#testmethods <- findMethods(computeValue)@names
testanalysis <- initializeanalysisclassobject(testobjectsfromparameters, testdataobject)
testreport <- initializeReportClass(testanalysis)
testrun <- preproviz(testdata)

# Test: DataClass object can be initialized
test_that("Object initialization", {expect_is(testdataobject, "DataClass")})

# Test: sub class object can be initialized as such or by giving ParameterClass and DataClass objects as arguments
test_that("Object initialization", {expect_is(testmissingvalueshare, "MissingValueShare")})
test_that("Object initialization", {expect_is(testobjectsfromparameters, "list")})

# Test: sub class object is valid
test_that("Object initialization", {expect_equal(testmissingvalueshareisvalid, FALSE)})

# Test: feature can be constructed and found as computeValue method
# test_that("Object initialization", {expect_equal("testfeature" %in% testmethods, TRUE)})

# Test: AnalysisClass object can be initialized
test_that("Object initialization", {expect_is(testanalysis, "AnalysisClass")})

# Test: ReportClass object can be initialized
test_that("Object initialization", {expect_is(testreport, "ReportClass")})

# Test: RunClass object can be initialized
test_that("Object initialization", {expect_is(testrun, "RunClass")})
