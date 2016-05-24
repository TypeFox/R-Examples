ivar <- function(x, ns=attr(x,"ns"), drop.levels=attr(x,"drop.levels")){
#
#   Rep a column vector into a matrix of correct size
#

if( is.null(ns) ) stop( paste("Number of occasions for", deparse(substitute(x)), "not specified"))

if( is.null(drop.levels) ) drop.levels <- 1

if(is.factor(x)){
    x <- as.factor(as.character(x))  # drop unused levels
    xm <- model.matrix( ~-1+x, data=x )   # no intercept, drop a level in F.cr.model.matrix if appropriate
    nan <- nrow(xm)  # nrow = n individuals, ncol = n contrast indicator vars
    nc <- ncol(xm)
    drop.levels <- drop.levels[ drop.levels>=1 & drop.levels <=nc ]  # drop.levels=0 means "do not drop"
    nd <- length(drop.levels)
    xfull <- NULL
    xlevs <- NULL
    for( j in 1:nc){
        if( all(j != drop.levels) ){
            xfull <- cbind(xfull, matrix( xm[,j], nrow=nan, ncol=ns ))
            xlevs <- c(xlevs, levels(x)[j])
        }            
    }
    attr(xfull,"levels") <- xlevs
    attr(xfull,"contr")<-attr(xm, "contrasts")[[1]]
} else {    
    xfull <- matrix( x, nrow=length(x), ncol=ns)
}
xfull

}
