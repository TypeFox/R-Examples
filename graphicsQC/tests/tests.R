# ------------------- Testing makeplots.R -------------------------
library(graphicsQC)
# General Tests:
# plotExpr takes arguments; expr, filetype, prefix, path, clear
test1 <- plotExpr(c("plot(1:10)", "plot(4:40)", "x<-3", "plot(2:23)"),
                  c("pdf", "ps", "png"), "testdir1/", "test", FALSE)
# There should now be a folder "testdir1" in getwd() containing
# pdf and ps files and test-log.xml, ie
list.files("testdir1")

test1a <- plotExpr(c("plot(1:10)", "warning(\"warn\")", "x<-3", "plot(2:23)"),
                   c("pdf", "ps", "png"), "test1a/", "test", FALSE)
test1b <- plotExpr(c("plot(1:10)", "x<-3", "plot(2:23)"),
                   c("pdf", "ps", "png"), "test1b/", "test", FALSE)
testA <- compare(test1a, test1b, erase = "none")

# No file result
blankTest <- plotExpr(c("x<-3", "x+2"), c("pdf", "ps", "png"),
                      "blankTest/", "test", FALSE)
blankTest2 <- plotExpr(c("x<-2", "x+3"), c("pdf", "png"),
                       "blankTest2/", "blankTest", FALSE)

# plotExpr() allows expressions as well as text
expressions <- expression({ plot(1:10); plot(4:40); x<-3; plot(2:23) })
test1Expr <- plotExpr(expressions, "png", "exprDir", "exprTest", FALSE)
list.files("exprDir")

# Repeating the same again should give an error as clear is False
try(plotExpr(c("plot(1:10)", "plot(4:40)", "x<-3", "plot(2:23)"),
                c("pdf","ps"), "testdir1", "test", FALSE))
# Error: files of intended filename already exist in ‘path’

# plotFunction takes same arguments but prefix defaults to fun
testFun <- plotFunction(c("plot", "lm"), c("pdf", "ps"), path="testdir2",
                      clear=FALSE)
# Folder with output is produced in testdir2, with plot-log.xml,
# lm-log.xml, and plot-lm-funLog.xml
list.files("testdir2")

# Test plotFile
plotFiletest <- plotFile(file.path("testFiles", c("Rfile.R", "Rfile2.R")),
                         "png", path="filedir1")
list.files("filedir1")

## --Insert tests for plotPackage when done

# Testing evalPlotCode with badExpression
badExpression <- c("x<-3", "plot(x:8)", "warning(\"firstWarning\")",
          "warning(\"secondWarning\")", "stop(\"end error\")",
          "warning(\"NOTCALLED\")")
lapply(c("pdf", "ps"), graphicsQC:::evalPlotCode, badExpression, "prefix-%d")

# ------------------- Testing compare.R -------------------------
# --- Testing reading of XML log files ---
test1Check <- graphicsQC:::readLog(file.path("testdir1",
                                             "test-log.xml"))
identical(test1, test1Check) # Should be True

testBlankCheck <- graphicsQC:::readLog(file.path("blankTest",
                                             "test-log.xml"))
identical(blankTest, testBlankCheck) # Should be True

testFunCheck <- graphicsQC:::readLog(file.path("testdir2",
                                               "plot-funLog.xml"))
identical(testFun, testFunCheck) # Should be True

plotFileCheck <- graphicsQC:::readLog(file.path("filedir1",
                                                "Rfile.R-fileLog.xml"))
identical(plotFiletest, plotFileCheck) # Should be True

testFunAutoCheck <- graphicsQC:::getQCResult("testdir2") # Testing auto-detect
identical(testFun, testFunAutoCheck)

## --Insert test for readPlotPackageLog when done

# --- Testing compare for plotExpr ---
blankComp <- compare(blankTest, blankTest2, erase="none")
# no comparisons/no files

blankCompCheck <- graphicsQC:::readLog(file.path("blankTest",
                                       "test+blankTest-compareExprLog.xml"))
identical(blankComp, blankCompCheck)
                                          
compare(test1, test1, erase="none") # all identical

# Same as test1 but missing the last plot and png
test1b <- plotExpr(c("plot(1:10)", "plot(4:40)", "plot(2:23)",
                     "warning(\"a\")", "warning(\"b\")", "warning(\"c\")"),
                   c("pdf", "ps"), "testdir1b/", "test", FALSE)
test1B <- plotExpr(c("plot(1:10)", "plot(4:40)", "warning(\"d\")",
                     "warning(\"a\")"), c("pdf", "png"), "testdir1Be", "test2")
test1BComp <- compare(test1b, test1B, erase="none") # Should all be identical
                                # but with unpaired files (test1b has an extra 
                                # plot and ps; test1B has png)
# Test round-trip for compare
test1BCompCheck <- graphicsQC:::readLog(file.path("testdir1b",
                                        "test+test2-compareExprLog.xml"))
identical(test1BComp, test1BCompCheck) # Should be True

# Second plot is different, third plot is still unpaired, test1B will be wiped
# with clear being TRUE
test1C <- plotExpr(c("plot(1:10)", "plot(4:41)"), c("pdf", "png"), 
                   "testdir1Be", "test2", TRUE)
compare(test1, test1C, erase="none") # This creates (and overwrites!) a
                                     # compareLog in the test dir, i.e.
list.files("testdir1") #test+test2-compareExprLog.xml is new

# --- Testing compare for plotFun ---
# All plots from plot are different to barplot
testFun1 <- plotFunction(c("plot", "plot"), c("pdf", "ps"), prefix=c("plot",
                         "plot2"), path="testFun", clear=FALSE)
testFun1ReadCheck <- graphicsQC:::readLog(file.path("testFun",
                                                    "plot-funLog.xml"))
identical(testFun1, testFun1ReadCheck)

testFun2 <- plotFunction(c("barplot", "plot"), c("pdf", "ps"), path="testFun2",
                         clear=FALSE)
testFun2ReadCheck <- graphicsQC:::readLog(file.path("testFun2",
                                                    "barplot-funLog.xml"))
identical(testFun2, testFun2ReadCheck)

funComparison <- compare(testFun1, testFun2, erase="none")
funComparisonReadCheck <- graphicsQC:::readLog(file.path("testFun",
                                               "plot-compareFunLog.xml"))
identical(funComparison, funComparisonReadCheck)

# --- Testing compare for plotFile ---
testFile3 <- plotFile(file.path("testFiles", c("Rfile.R", "Rfile3.R")),
                      c("pdf", "png"), prefix=c("file1pref", "file2pref"),
                      path="testFile3")
testFile4 <- plotFile(file.path("testFiles", c("Rfile2.R", "Rfile4.R")),
                      c("pdf", "png", "ps"), path="testFile4", clear = FALSE)
testFile5 <- plotFile(file.path("testFiles", c("Rfile.R", "Rfile2.R")),
                      c("pdf", "png", "ps"), path="testFile5")
fileComparison <- compare(testFile3, testFile4, erase="none")
fileComparisonReadCheck <- graphicsQC:::readLog(file.path("testFile3",
                                               "file1pref-compareFileLog.xml"))
identical(fileComparison, fileComparisonReadCheck)

# Note: Diff plots are being stored in the test directory, so
list.files("testFun")

testFun3 <- plotFunction(c("plot", "barplot"), c("pdf", "png"),
                         path="testFun3", clear=FALSE)
list.files("testFun3")

# ------------------- Testing report.R -------------------------
if (require("Sxslt")) {
    report1 <- writeReport(test1, browse=FALSE)
    report2 <- writeReport(test1BComp, browse=FALSE)
    report3 <- writeReport(blankComp, browse=FALSE)
    report4 <- writeReport(funComparison, browse=FALSE)
    report5 <- writeReport(fileComparison, browse=FALSE)
    testFun4 <- plotFunction(c("barplot", "barplot"), c("pdf", "ps"),
                             prefix=c("bplot", "bplot2"),
                             path="testFun4", clear=FALSE)
    report6 <- writeReport(compare(testFun4, testFun1, "report6",
                           erase="none"), browse=FALSE)
    writeReport(testA, browse=FALSE)
}

## Something to note: Doing something like
# y <- 1:10
# x <- rnorm(10)
# myLm <- lm(y~x)
# plotExpr("plot(myLm)", "png", "testlm", "testlm", FALSE)
# Will require the user to keep doing "Hit <Return> to see next plot:"
# if not in batch mode
# FULL fix requires change to source code of plot.lm
# WORKAROUND is just use BATCH mode !



