## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(stubthat)

## ------------------------------------------------------------------------
jedi_or_sith <- function(x) return('No one')
jedi_or_sith_stub <- stub(jedi_or_sith)

## ------------------------------------------------------------------------
jedi_or_sith_stub$withArgs(x = 'Luke')$returns('Jedi')

## ------------------------------------------------------------------------
jedi_or_sith('Luke')
jedi_or_sith_stub$f('Luke')

## ------------------------------------------------------------------------
library('httr')

check_api_endpoint_status <- function(url) {
  response <- GET(url)
  response_status <- status_code(response)
  ifelse(response_status == 200, 'up', 'down')
}

## ------------------------------------------------------------------------
stub_of_get <- stub(GET)
stub_of_get$withArgs(url = 'good url')$returns('good response')
stub_of_get$withArgs(url = 'bad url')$returns('bad response')

stub_of_status_code <- stub(status_code)
stub_of_status_code$withArgs(x = 'good response')$returns(200)
stub_of_status_code$withArgs(x = 'bad response')$returns(400)

library('testthat')
with_mock(GET = stub_of_get$f, status_code = stub_of_status_code$f,
          expect_equal(check_api_endpoint_status('good url'), 'up'))

with_mock(GET = stub_of_get$f, status_code = stub_of_status_code$f,
          expect_equal(check_api_endpoint_status('bad url'), 'down'))

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$expects(a = 2)

stub_of_identify$f(2)
stub_of_identify$f(3)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$strictlyExpects(a = 2)
stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
stub_of_identify$strictlyExpects(a = 2, b = 1)
stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(3)$expects(a = 2)

stub_of_identify$f(100)
stub_of_identify$f(100)
stub_of_identify$f(100)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(3)$strictlyExpects(a = 2, b = 2)

stub_of_identify$f(2)
stub_of_identify$f(2)
stub_of_identify$f(2)

## ------------------------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$returns(0)

stub_of_identify$f(2)

## ------------------------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(2)$returns(0)

stub_of_identify$f(2)
stub_of_identify$f(2)

## ------------------------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$withArgs(a = 2)$returns(0)

stub_of_identify$f(1)
stub_of_identify$f(2)

## ------------------------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$withExactArgs(a = 2)$returns(0) # won't work because value for b is not defined
stub_of_identify$withExactArgs(a = 2, b = 1)$returns(1)

stub_of_identify$f(1)
stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$throws('some err msg')

stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(2)$throws('some err msg')

stub_of_identify$f(0)
stub_of_identify$f(0)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$withArgs(a = 2)$throws('some err msg')

stub_of_identify$f(1)
stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$withExactArgs(a = 2)$throws('good') # won't work because value for b is not defined
stub_of_identify$withExactArgs(a = 2, b = 1)$throws('nice')

stub_of_identify$f(1)
stub_of_identify$f(2)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

ans <- stub_of_identify$f(3)
ans <- stub_of_identify$f(3)
stub_of_identify$calledTimes()
ans <- stub_of_identify$f(3)
stub_of_identify$calledTimes()

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(1)$expects(a = 1)$returns('good')
stub_of_identify$onCall(3)$expects(a = 3)$returns('nice')

stub_of_identify$f(3)
stub_of_identify$f(3)
stub_of_identify$f(3)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(1)$expects(a = 1)
stub_of_identify$onCall(1)$returns('good')
stub_of_identify$onCall(3)$returns('nice')
stub_of_identify$onCall(3)$expects(a = 3)

stub_of_identify$f(3)
stub_of_identify$f(3)
stub_of_identify$f(3)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(1)$strictlyExpects(a = 3)$returns('good')
stub_of_identify$onCall(3)$strictlyExpects(a = 3, b = 1)$returns('nice')

stub_of_identify$f(3)
stub_of_identify$f(3)
stub_of_identify$f(3)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(1)$expects(a = 1)$throws('good')
stub_of_identify$onCall(3)$expects(a = 3)$throws('nice')

stub_of_identify$f(3)
stub_of_identify$f(3)
stub_of_identify$f(3)

## ----error=TRUE----------------------------------------------------------
sum <- function(a, b = 1) return(a + b)
stub_of_identify <- stub(sum)

stub_of_identify$onCall(1)$strictlyExpects(a = 3)$throws('good')
stub_of_identify$onCall(3)$strictlyExpects(a = 3, b = 1)$throws('nice')

stub_of_identify$f(3)
stub_of_identify$f(3)
stub_of_identify$f(3)

