library(RUnit)

test_that(
  "checkTrue gets converted to expect_true",
  {
    test_truth <- function()
    {
      x <- all(runif(10) > 0)
      checkTrue(x)
    } 
    expected <- quote(
      test_that(
        "test_truth",
        {
          x <- all(runif(10) > 0)
          expect_true(x)
        }
      )    
    )
    actual <- convert_test(test_truth)
    expect_equal(actual, expected)    
  }
)

test_that(
  "not checkTrue gets converted to expect_false",
  {
    test_falsity <- function()
    {
      x <- any(runif(10) < 0)
      checkTrue(!x)
    } 
    expected <- quote(
      test_that(
        "test_falsity",
        {
          x <- any(runif(10) < 0)
          expect_false(x)
        }
      )    
    )
    actual <- convert_test(test_falsity)
    expect_equal(actual, expected)    
  }
)

test_that(
  "checkEquals gets converted to expect_equals",
  {
    test_equality <- function()
    {
      x <- sqrt(1:5)
      expected <- c(1, 4, 9, 16, 25)
      checkEquals(expected, x ^ 4)
    } 
    expected <- quote(
      test_that(
        "test_equality",
        {
          x <- sqrt(1:5)
          expected <- c(1, 4, 9, 16, 25)
          expect_equal(x ^ 4, expected)
        }
      )    
    )
    actual <- convert_test(test_equality)
    expect_equal(actual, expected)    
  }
)

test_that(
  "checkEqualsNumeric gets converted to expect_equals",
  {
    test_equality <- function()
    {
      x <- sqrt(1:5)
      expected <- c(1, 4, 9, 16, 25)
      checkEqualsNumeric(expected, x ^ 4)
    } 
    expected <- quote(
      test_that(
        "test_equality",
        {
          x <- sqrt(1:5)
          expected <- c(1, 4, 9, 16, 25)
          expect_equal(x ^ 4, expected)
        }
      )    
    )
  actual <- convert_test(test_equality)
  expect_equal(actual, expected)    
  }
)

test_that(
  "checkIdentical gets converted to expect_identical",
  {
    test_identicality <- function()
    {
      x <- paste(letters[1:10], collapse = " ")
      expected <- "a b c d e f g h i j"
      checkIdentical(expected, x)
    } 
    expected <- quote(
      test_that(
        "test_identicality",
        {
          x <- paste(letters[1:10], collapse = " ")
          expected <- "a b c d e f g h i j"
          expect_identical(x, expected)
        }
      )    
    )
    actual <- convert_test(test_identicality)
    expect_equal(actual, expected)    
  }
)

test_that(
  "checkException gets converted to expect_error",
  {
    test_error <- function()
    {
      checkException("1" + "2")
    }
    expected <- quote(
      test_that(
        "test_error",
        {
          expect_error("1" + "2")
        }
      )
    )
    actual <- convert_test(test_error)
    expect_equal(actual, expected)    
  }
)

test_that(
  "simple, unbraced check* functions are converted to expect_* functions.",
  {
    # Notice that the arguments are (deliberately) swapped in checkEquals vs. 
    # expect_equal.
    runit_test <- function()
    {
      x <- 1:10
      checkTrue(all(x > 0))
      checkTrue(!any(x < 0))
      checkIdentical(1:10, x)
      checkEquals(seq.int(2, 20, 2) / 2, x)
      checkEqualsNumeric(seq.int(2, 20, 2) / 2, x)
      checkException(stop("!!!"))
    }   
    expected <- quote(
      test_that(
        "runit_test",
        {
          x <- 1:10
          expect_true(all(x > 0))
          expect_false(any(x < 0))
          expect_identical(x, 1:10)
          expect_equal(x, seq.int(2, 20, 2) / 2)
          expect_equal(x, seq.int(2, 20, 2) / 2)
          expect_error(stop("!!!"))
        }
      )  
    )
    expect_equal(convert_test(runit_test, "runit_test"), expected)
  }
)

test_that(
  "checks inside braces are converted.",
  {
    test_braces <- function()
    {
      {
        {
          x <- 2
          checkTrue(x > 1)
        }
      }
    }   
    expected <- quote(
      test_that(
        "test_braces",
        {
          {
            {
              x <- 2
              expect_true(x > 1)
            }
          }
        }
      )
    )
    actual <- convert_test(test_braces)
    expect_equal(actual, expected)
  }
)

test_that(
  "checks inside for loops (no braces) are converted.",
  {
    test_for_loops_no_braces <- function()
    {
      for(i in 1:3) checkTrue(i > 0)
    }   
    expected <- quote(
      test_that(
        "test_for_loops_no_braces",
        {
          for(i in 1:3) expect_true(i > 0)
        }
      )
    )
    actual <- convert_test(test_for_loops_no_braces)
    expect_equal(actual, expected)
  }
)

test_that(
  "checks inside for loops (with braces) are converted.",
  {
    test_for_loops_with_braces <- function()
    {
      for(i in 1:3)
      {
        x <- i
        checkTrue(x > 0)
      }
    }   
    expected <- quote(
      test_that(
        "test_for_loops_with_braces",
        {
          for(i in 1:3)
          {
            x <- i
            expect_true(x > 0)
          }
        }
      )
    )
    actual <- convert_test(test_for_loops_with_braces)
    expect_equal(actual, expected)
  }
)

test_that(
  "checks inside switch statements are converted.",
  {
    test_switch <- function()
    {
      switch(
        "b", 
        a = mean(1:5),
        b = checkEquals(0, 0)
      )
    }   
    expected <- quote(
      test_that(
        "test_switch",
        {
          switch(
            "b", 
            a = mean(1:5),
            b = expect_equal(0, 0)
          )
        }
      )
    )
    actual <- convert_test(test_switch)
    expect_equal(actual, expected)
  }
)

test_that(
  "checks inside nested switch statements are converted.",
  {
    test_nested_switch <- function()
    {
      x <- "richie"
      switch(
        x,
        richie = {
          switch(
            x,
            richie = {
              checkException(stop("!!!"))
            }
          )
        },
        cotton = {
          checkIdentical("cotton", x)
        }
      )
    }   
    expected <- quote(
      test_that(
        "test_nested_switch",
        {
          x <- "richie"
          switch(
            x,
            richie = {
              switch(
                x,
                richie = {
                  expect_error(stop("!!!"))  
                }
              )
            },
            cotton = {
              expect_identical(x, "cotton")
            }
          )
        }
      )
    )
    actual <- convert_test(test_nested_switch)
    expect_equal(actual, expected)
  }
)

test_that(
  "checks inside (possibly nested) loops or if blocks are converted.",
  {
    test_nested_loops <- function()
    {
      x <- 6:10
      for(i in 1:5)
      {
        if(i %% 2 == 0)
        {
          checkTrue(all(x > i), msg = "i divisible by 2") 
          if(i == 4)
          {
            checkIdentical(4, i, msg = "i = 4")
          } else
          {
            while(i > 0) 
            {
              checkIdentical(2, i, msg = "i = 2")
            }
            repeat
            {
              checkException(stop("!!!"))
              break
            }
          }
        }               
      }
    }
    expected <- quote(
      test_that(
        "test_nested_loops",
        {
          x <- 6:10
          for(i in 1:5)
          {
            if(i %% 2 == 0)
            {
              expect_true(all(x > i), info = "i divisible by 2") 
              if(i == 4)
              {
                expect_identical(i, 4, info = "i = 4")
              } else
              {
                while(i > 0) 
                {
                  expect_identical(i, 2, info = "i = 2")
                }
                repeat
                {
                  expect_error(stop("!!!"))
                  break
                }
              }
            }               
          }
        }
      )  
    )
    actual <- convert_test(test_nested_loops)
    expect_equal(actual, expected) 
  }
)

test_that(
  "failing tests still return a testthat test.",
  {
    # This is will not be the case if the RUnit test is in a file (due to the way
    # that sys.source works).
    failing_test <- function()
    {
      checkEquals(1, 2)
    }
    expected <- quote(
      test_that(
        "failing_test",
        {
          expect_equal(2, 1)
        }
      )  
    )
    actual <- convert_test(failing_test)
    expect_equal(actual, expected) 
  }
)


