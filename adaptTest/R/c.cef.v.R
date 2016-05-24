`c.cef.v` <-
function (c)        {if (is.na(c)) {result <- NA} else {if(!(le(0,c)))                            {result <- NA} else {if(eq(c,0)) {result <- function (t) 0+t-t} else {if(eq(c,Inf)) {result <- function (t) 1+t-t} else {result <- function (t) (1-t^c)^(1/c)}}}}; return(result)}

