scacum <- structure(function#Scaled-cummulative sums
### Cummulative sums of a vector are scaled on a constant value. 
##details<< Cummulative values of the numeric vector are computed with \code{\link{cumsum}} 
(
    x, ##<<\code{numeric} vector with vector names being ordered on time
    y   =   NA,   ##<<\code{NA} or \code{numeric} constant to scale the
                  ##vector.  If \code{NA} then the cummulative sums
                  ##are not scaled.
    z = NA ##<<\code{NA} or \code{numeric} constant in range of
           ##the vector names. If NA then maximun value in such a
           ##range is used.
    
) {
    
    csn. <- FALSE
    if(is.data.frame(x)){
        csnu <- colclass(x,TRUE)[['num']]
        csn <- c(colclass(x,TRUE)[c('tmp','fac')],
                 recursive = TRUE)
        csn. <- length(csn)!=0
        csn.. <- csn[!csn%in%'csx']
        cd <- x
        x <- x[,'x']
        names(x) <- cd[,'year']}
        
    if(is.null(names(x)))
        stop('NULL year names of x',
             call. = FALSE)
    xcum <- cumsum(x)
    if(is.na(z))
        z <- max(as.numeric(names(x)))
    inc <- 0
    if(!is.na(y))
        inc <- y - xcum[as.character(z)]
    csx <- xcum + inc
    if(any(csx < 0,na.rm = TRUE))
        csx <- xcum
    xd <- data.frame(x,csx)
    
    if(csn.&& length(csnu) > 1){
        xd <- cd[,csnu]
        xd[,'csx'] <- csx }
    if(csn.)
        xd <- cbind(xd,cd[,csn..])
    
    return(xd)
### data frame with the original vector, and its scaled-cummulative sums.
} , ex=function() {
    x <- c(0.79,0.32,0.53,0.43,0.18)
    names(x) <- 1948:1952
    y <- 4
    z <- 1951
    scacum(x,y,z)
    
    ##If y = NA then cummulative values are scaled arround
    ##max(cumsum(x)):
    max(cumsum(x))
    scacum(x,NA)
})
