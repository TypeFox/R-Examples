tvar <- function(x, nan=attr(x,"nan"), drop.levels=attr(x,"drop.levels")){
#
#   Rep a row vector into a matrix of correct size
#

if( is.null(nan) ) stop( paste("Number of individuals for", deparse(substitute(x)), "not specified"))

if( is.null(drop.levels) ) drop.levels <- 1

if(is.factor(x)){
    x <- as.factor(as.character(x))  # drop unused levels
    xm <- model.matrix( ~-1+x, data=x )  
    ns <- nrow(xm)  # nrow = n occasions, ncol = n contrast indicator vars
    nc <- ncol(xm)
    drop.levels <- drop.levels[ drop.levels>=1 & drop.levels <=nc ]  # drop.levels=0 means "do not drop"
    nd <- length(drop.levels)
    xfull <- NULL
    xlevs <- NULL
    for( j in 1:nc){
        if( all(j != drop.levels) ){
            xfull <- cbind( xfull, matrix( xm[,j], nrow=nan, ncol=ns, byrow=TRUE ))
            xlevs <- c(xlevs, levels(x)[j])
        }
    }
    attr(xfull,"levels") <- xlevs
    attr(xfull,"contr") <- attr(xm, "contrasts")[[1]]
} else {    
    xfull <- matrix( x, nrow=nan, ncol=length(x), byrow=TRUE)
}
xfull

}
