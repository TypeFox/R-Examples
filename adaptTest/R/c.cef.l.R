`c.cef.l` <-
function (c)        {if (is.na(c)) {result <- NA} else {                                                               if(eq(c,Inf)) {result <- function (t) 0+t-t} else {if(eq(c,-Inf)) {result <- function (t) 1+t-t} else {result <- function (t) pnorm(-sqrt(2)*c-qnorm(t))}}}; return(result)}

