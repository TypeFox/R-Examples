context("Binary")

test_that("Rational example works", {
  suppressWarnings({
    Rational <- defineClass("Rational", contains = c("Show", "Binary"), {
      
      numer <- 0
      denom <- 1
      .g <- 1
      
      .gcd <- function(a, b) if(b == 0) a else Recall(b, a %% b)
      
      init <- function(numer, denom) {
        .self$.g <- .gcd(numer, denom)
        .self$numer <- numer / .g
        .self$denom <- denom / .g
      }
      
      show <- function() {
        cat(paste0(.self$numer, "/", .self$denom, "\n"))
      }
      
      ".+" <- function(that) {
        Rational(numer = numer * that$denom + that$numer * denom,
                 denom = denom * that$denom)
      }
      
      neg <- function() {
        Rational(numer = -.self$numer,
                 denom = .self$denom)
      }
      
      ".-" <- function(that) {
        .self + that$neg()
      }
      
    })
    
    rational <- Rational(2, 3)
    
    expect_equal(rational$numer, 2)
    expect_equal(rational$denom, 3)
    
    rational <- rational + rational
    expect_equal(rational$numer, 4)
    expect_equal(rational$denom, 3)
    
    rational <- rational$neg()
    expect_equal(rational$numer, -4)
    expect_equal(rational$denom, 3)
    
    rational <- rational - rational
    expect_equal(rational$numer, 0)
    expect_equal(rational$denom, 1)
    
    x <- Rational(numer = 1, denom = 3)
    y <- Rational(numer = 5, denom = 7)
    z <- Rational(numer = 3, denom = 2)
    
    rational <- x - y - z
    expect_equal(rational$numer, -79)
    expect_equal(rational$denom, 42)
  })
  
})
