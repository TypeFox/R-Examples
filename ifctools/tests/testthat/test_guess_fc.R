context("guess_fc tests")

test_that("Regular fiscal codes", {
            ## using fictious data
            expect_identical(guess_fc("Rossi",
                                      "Mario",
                                      as.Date("1960-01-01"),
                                      FALSE,
                                      "F205"
                                      ),
                             "RSSMRA60A01F205T")
            expect_identical(guess_fc("Bianchi",
                                      "Giovanna",
                                      as.Date("1970-01-01"),
                                      TRUE,
                                      "H501"
                                      ),
                             "BNCGNN70A41H501V")
          })

