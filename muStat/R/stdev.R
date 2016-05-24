`stdev` <-
function(x,na.rm,unbiased)
    sd(x,na.rm=TRUE)*if (unbiased) 1 else 1-1/length(x)

