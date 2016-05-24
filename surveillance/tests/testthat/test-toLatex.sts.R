context("toLatex-method for the sts-class")

data("ha.sts")
data("salmonella.agona")

test_that("toLatex accepts basic input and returns Latex", {  
  control <- list(
    noPeriods=10,populationBool=FALSE,
    fitFun="algo.farrington.fitGLM.flexible",
    b=4,w=3,weightsThreshold=2.58,
    pastWeeksNotIncluded=26,
    pThresholdTrend=1,trend=TRUE,
    thresholdMethod="new",alpha=0.01
  )
  result <- ha.sts
  result@alarm[,7]  <- TRUE
  result@upperbound[,7]  <- 1
  laTex <- toLatex(result, subset=(280:290), table.placement="h", size = "scriptsize",
                       sanitize.text.function = identity,
                       NA.string = "-",include.rownames=FALSE)
    
  laTex3 <- toLatex(result, subset=(280:290),
                    alarmPrefix = "aaaa",
                    alarmSuffix = "bbbb", table.placement="h", size = "scriptsize",
                   sanitize.text.function = identity,
                   NA.string = "-",include.rownames=FALSE)

  expect_true(grepl("aaaa", paste(as.character(laTex3), collapse = ' ')))
  expect_true(grepl("bbbb", paste(as.character(laTex3), collapse = ' ')))
  expect_is(laTex, "Latex")
  expect_is(laTex3, "Latex")
})

test_that("caption is incorporated", {
  testCaption <- "Please print my caption"
  latex <- toLatex(ha.sts, caption = testCaption)
  expect_true(grepl(testCaption, paste(as.character(latex), collapse = ' ')))
})

test_that("label is incorporated", {
  testLabel <- "Please print my label"
  latex <- toLatex(ha.sts, label = testLabel)
  expect_true(grepl(testLabel, paste(as.character(latex), collapse = ' ')))
})

test_that("ubColumnLabel is incorporated", {
  testUBLabel <- "Upperbound"
  latex <- toLatex(ha.sts, ubColumnLabel = testUBLabel)
  expect_true(grepl(testUBLabel, paste(as.character(latex), collapse = ' ')))
})

test_that("one can override the default table column labels", {
  columnLabels <- c("Jahr", "Woche", "chwi1", "UB",
                    "frkr2", "UB", "lich3", "UB", "mahe4", "UB", "mitt5", "UB", 
                    "neuk6", "UB", "pank7", "UB", "rein8", "UB", "span9", "UB", 
                    "zehl10", "UB", "scho11", "UB", "trko12", "UB")
  latex <- toLatex(ha.sts, columnLabels = columnLabels)
  expect_true(all(
    sapply(columnLabels, 
               function(l) grepl(l, paste(as.character(latex), collapse = ' '))
           , USE.NAMES = FALSE) 
    ))
})

test_that("toLatex works with output from farringtonFlexible()", {
  # Create the corresponding sts object from the old disProg object
  salm <- disProg2sts(salmonella.agona)
  # Farrington with old options
  control1 <- list(range=(260:312),
                   noPeriods=1,populationOffset=FALSE,
                   fitFun="algo.farrington.fitGLM.flexible",
                   b=4,w=3,weightsThreshold=1,
                   pastWeeksNotIncluded=3,
                   pThresholdTrend=0.05,trend=TRUE,
                   thresholdMethod="delta",alpha=0.1)
  salm1 <- farringtonFlexible(salm,control=control1)  
  expect_is(toLatex(salm1), "Latex")
})

test_that("toLatex only accepts a single sts object", {
  expect_error(toLatex(list(ha.sts, ha.sts)))
})

test_that("toLatex stops if 'subset' is not applicable", {
  expect_error(toLatex(ha.sts, subset=(-5:290)))
  expect_error(toLatex(ha.sts, subset=(1:10000)))
  expect_error(toLatex(ha.sts, subset=(10000:100000)))
})
