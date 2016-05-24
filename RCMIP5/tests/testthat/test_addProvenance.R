# Testing code for the RCMIP5 scripts in 'addProvenance.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("addProvenance")

test_that("addProvenance handles bad input", {
    expect_error(addProvenance())
    expect_error(addProvenance(x=3))
    expect_error(addProvenance(cmip5data(), msg=3))    
})

test_that("addProvenance initializes", {
    d <- addProvenance(cmip5data(1), msg="test")
    expect_is(d, "cmip5data")
    expect_is(d$provenance, "data.frame")
})

test_that("addProvenance adds messages", {
    d <- cmip5data(1)
    expect_equal(nrow(d$provenance), 2)            # One line for software specs, one for message
    d <- addProvenance(d, "test23")
    expect_equal(nrow(d$provenance), 3)           # One line for software specs, one for message
    expect_true(any(grepl("test23", d$provenance[3,])))   # 'test23' should appear in final line
})

test_that("addProvenance merges provenances", {
    d1 <- cmip5data(1)
    d2 <- cmip5data(2)
    d3 <- addProvenance(d1, d2)
    
    expect_equal(nrow(d1$provenance)+nrow(d2$provenance), nrow(d3$provenance))
})
