contr.FrF2 <- function (n, contrasts=TRUE) 
{
    ## the contrasts option does not do anything but is needed for model.matrix
    ## to work on objects of this type
    if (!contrasts) stop("contr.FrF2 not defined for contrasts=FALSE")

    ## CAUTION: for more than 4 levels, levels need to be in correct order
    ## for the FrF2 structure to hold!
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("invalid choice for n in contr.blocks")
    }
    else levels <- n
    lenglev <- length(levels)
    if (!2^round(log2(lenglev))==lenglev) 
        stop("contr.FrF2 requires that the number of levels is a power of 2.")

    ## definition of contrast matrix
       if (lenglev==2) destxt <- "matrix(c(-1,1),ncol=1)"
       else {
       destxt <- "expand.grid(c(-1,1)"
       for (i in 2:round(log2(lenglev))) 
                destxt <- paste(destxt,",c(-1,1)",sep="")
       destxt <- paste("as.matrix(",destxt,"))",sep="")
       }
       cont <- eval(parse(text=destxt))
       cont <- sapply(Yates[1:(lenglev-1)], function(obj) (apply(cont[,obj,drop=FALSE],1,prod)))
       rownames(cont) <- levels
       colnames(cont) <- NULL
    cont
}
