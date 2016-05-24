`c.cef.h` <-
function (c)        {if (is.na(c)) {result <- NA} else {if(!(le(0,c) && le(c,1)))                 {result <- NA} else {result <- function (t) c+t-t}}; return(result)}

