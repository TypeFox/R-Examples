`c.cef.b` <-
function (c)        {if (is.na(c)) {result <- NA} else {if(!(le(0,c) && le(c,1)))                 {result <- NA} else {if(eq(c,0)) {result <- function (t) 0+t-t} else {if(eq(c,1)) {result <- function (t) 1+t-t} else {result <- function (t) pmin(1, c/t)}}}}; return(result)}

