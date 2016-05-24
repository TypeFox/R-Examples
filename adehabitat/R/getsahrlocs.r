"getsahrlocs" <- function(x, what=c("sa", "hr", "locs"))
{
    ## Verifications
    if (!inherits(x, "sahrlocs"))
        stop("x should be of class sahrlocs")
    what<-match.arg(what)

    ## Core of the function
    sahr<-x
    rm(x)
    if (!inherits(sahr, "sahrlocs")) stop("non convenient data type")
    if (is.na(match(what, c("sa", "hr", "locs"))))
        stop("what should be either \"sa\", \"hr\", or \"locs\"")
    output<-sahr[[what]]


    ## Output
    attr(output, "nrow")<-attr(sahr, "nrow")
    attr(output, "ncol")<-attr(sahr, "ncol")
    attr(output, "xll")<-attr(sahr, "xll")
    attr(output, "yll")<-attr(sahr, "yll")
    attr(output, "cellsize")<-attr(sahr, "cellsize")
    class(output)<-c("kasc", "data.frame")

    return(output)
  }

