context("Message reporting and output levels")

test_that("message reporting works", {
    library(ore)
    
    options(reportrOutputLevel=NULL)
    options(reportrStderrLevel=OL$Fatal)
    
    expect_output(getOutputLevel(), "Output level is not set", fixed=TRUE)
    options(reportrOutputLevel=OL$Debug)
    expect_equivalent(getOutputLevel(), OL$Debug)
    setOutputLevel(OL$Info)
    expect_equivalent(getOutputLevel(), OL$Info)
    
    setOutputLevel(OL$Warning)
    expect_silent(report(OL$Info,"Test message"))
    
    setOutputLevel(OL$Info)
    expect_output(report(OL$Info,"Test message"), "Test message", fixed=TRUE)
    expect_output(report(Info,"Test message"), "Test message", fixed=TRUE)
    expect_output(report("Info","Test message"), "Test message", fixed=TRUE)

    flag(OL$Warning, "Test warning")
    flag(OL$Warning, "Test warning")
    expect_output(reportFlags(), "[x2]", fixed=TRUE)
    
    expect_warning(sqrt(-1), "NaNs produced", fixed=TRUE)
    expect_output(withReportrHandlers(sqrt(-1)), "NaNs produced", fixed=TRUE)
    
    f <- function() message("Howdy")
    expect_output(withReportrHandlers(f()), "* INFO: Howdy", fixed=TRUE)
    
    setOutputLevel(OL$Debug)
    options(reportrStackTraceLevel=OL$Info)
    expect_output(flag(Info,"Converted to report"), "Converted to report", fixed=TRUE)
    expect_output(withReportrHandlers(f()), "* f()", fixed=TRUE)
    
    setOutputLevel(OL$Error)
    expect_null(ask("This question will not be answered"))
})
