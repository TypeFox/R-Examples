rtimes <- structure(function#Relative times
### Relative times from 1 to the number of observed years are computed. Duplicated times can be sinchronized by replacing unique values with NA.   
(
    x, ##<< \code{numeric} vector with names representing the
       ##formation years, or multilevel data frame containing a column
       ##of years.
    only.dup = TRUE ##<< \code{logical}. Extract only duplicated
                    ##times.  If TRUE then unique times are replaced
                    ##with NA. If all computed times are unique then
                    ##this argument is ignored.
) {
    csn. <- FALSE
    if(is.data.frame(x)){
        csnu <- colclass(x,T)[['num']]
        csn <- c(colclass(x,T)[
                             c('tmp','fac')],recursive = T)
        csn. <- length(csn)!=0
        csn.. <- csn[!csn%in%'time']
        cd <- x
        x <- x[,'x']
        names(x) <- cd[,'year']}
    
    n <- as.numeric(names(x))
    time <- abs(min(n) - n - 1)
    da <- data.frame(x,time)
    ## return(csn.&& length(csnu) > 1)
    if(csn.&& length(csnu) > 1){
        ## csn.. <- csn[!csn%in%'time'] #
        da <- cd[,csnu]
        da[,'time'] <- time }
    rownames(da) <- 1:nrow(da)
    dp <- duplicated(da[,'time'])
    uni <- with(da,!time%in%da[,'time'][dp])
    if(only.dup&any(dp))
        da[uni,'time'] <- NA
    if(csn.)
        da <- cbind(da,cd[,csn..])
    return(da)
### \code{data.frame} object with the initial vector and its formation times.
} , ex=function(){
    ## row names of a vector
    fy <- function(y,span){(y - span):y}
    x <- c(fy(2005,5),fy(2007,10)) 
    ## (not run) Simulating the vector
    r <- abs(rnorm(length(x)))
    names(r) <- x
    ## (not run) computing relative times:
    rtimes(r,only.dup = TRUE)        
    ## (not run) Extracting only duplicated times:
    na.omit(rtimes(r,only.dup = TRUE))
})
