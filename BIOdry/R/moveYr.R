moveYr <- structure(function#Seasonal years
### Level of years is processed for each of the years can begin in a month different than 'January'.
##details<<\code{character} months as defined in \code{\link{month.abb}} or \code{\link{month.name}}.  
(
    cd, ##<<\code{data.frame} object with factor-level columns of
        ##months and years, or \code{numeric} vector of repeated years
        ##with vector names being the months.
    ini.mnt = 'Oct' ##<<\code{character} of initial month of the year
                    ##(see details), or \code{numeric} value from 1 to
                    ##12. Default 'Oct' produces years to begin in
                    ##October.
) {
    
    chn <- function(tmp){
        tmp. <- as.character(tmp)
        nm <- 1:12
        fm <- 'month.abb'
        mna <- tmp%in%month.abb[nm]
        if(!all(mna))
            fm <- 'month.name'
        names(nm) <- get(fm)[nm]
        if(!is.numeric(tmp)){
            tmp <- nm[tmp][tmp.]}
        names(tmp) <- get(fm)[tmp]
        return(tmp)}
    
    isdf <- is.data.frame(cd)
    if(isdf){
        ny <- cd
        cd <- cd[,'year']
        names(cd) <- chn(ny[,'month'])}
    if(!isdf)
        names(cd) <- chn(names(cd))
    ini.mnt <- chn(ini.mnt)
    mn <- 1:12
    ncd <- ifelse(
        mn[as.numeric(names(cd))] >= mn[ini.mnt],
           ifelse(
               mn[ini.mnt] > mn[5],
               cd + 1, cd),ifelse(
                               mn[ini.mnt] <= mn[5],
                               cd - 1, cd)) 
    if(isdf){
        ny[,'year'] <- ncd
        ny[,'month'] <- as.numeric(names(cd))}    
    if(!isdf){
        ny <-  ncd
        names(ny) <- month.abb[as.numeric(names(cd))]}
    return(ny)
### \code{data.frame} object with the months being \code{numeric}
### values and the years beguinning in \code{ini.mnt}
} , ex=function() {
    ## Climatic records of monthly precipitations and temperatures
    data(PTclim05,envir = environment())
    
    ## Making year 1955 in plot 'P16106'to begin on 'April'
    cl1 <- splitFrame(PTclim05,'year')[['P16106.1955']]
    cl2 <- moveYr(cl1,ini.mnt = 4)
    head(cl2)
    
    ## a simple vector of years
    yr <- rep(2005,12)
    names(yr) <- month.abb[1:12]
    moveYr(yr)
    
})
