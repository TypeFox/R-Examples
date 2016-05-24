ringApply <- structure(function#Multilevel apply
### Wrapper of \code{\link{Map}} to apply functions on multilevel data frames and preserve factor-level structure in the outputs.
##details<< Other functions such as \code{\link{rtimes}}, \code{\link{scacum}}, \code{\link{amod}}, or \code{\link{wlai}} can be implemented.  Function arguments should be formulated as suggested in \code{\link{mapply}}, with no vectorized arguments being stored in a \code{MoreArgs} list.
##references<< Lara, W., F. Bravo, D. Maguire. 2013. Modeling patterns between drought and tree biomass growth from dendrochronological data: A multilevel approach. Agric. For. Meteorol., 178-179:140-151.                       
(
    rd, ##<<\code{data.frame} object with factor-level columns.
    lv = 1, ##<< {numeric} position, or {character} name, of a
    ##factor-level column in the data frame.
    fn = 'scacum', ##<< \code{character} name of the function to be
    ##evaluated (see details). Default 'scacum'
    ##computes scaled-cummulative radii.
    ... ##<< Further arguments in \code{fn} (see details)
) {
    
    levs <- colclass(rd,TRUE)[['fac']]
    if(is.numeric(lv)) lv <- levs[lv]
    
    cl1 <- splitFrame(rd,lv)
    
    fam <- function(x,...){
        do.call(fn,list(x,...))}
    
    cl2 <- Map(function(x,...)fam(x,...), cl1,...)
    nord <- names(cl2)[order(names(cl2))]
    cl3 <- cl2[nord]
    cl4 <- do.call(rbind,cl3)
    rownames(cl4) <- NULL
    return(cl4)
### \code{data.frame} object preserving initial factor-level columns.
} , ex=function() {

    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    ## Radial increments measured on 2003:
    data(Pradii03,envir = environment())    
    ## Monthly precipitation sums and average temperatures:
    data(PTclim05,envir = environment())
    
    ##synchronized times by tree with 'rtimes':
    dfm1 <- ringApply(Prings05,lv = 2,fn = 'rtimes')
    str(dfm1)
    ##some synchronized times from time 1 to time 9:
    subset(dfm1,time%in%c(1:9,NA))
    
    ## Cummulative radial increments by sample:
    dfm2 <- ringApply(dfm1,lv = 'sample',y = Pradii03,fn = 'scacum')
    str(dfm2)    
    ##Allometric modeling by sample:
    dfm3 <- ringApply(dfm2,lv = 'sample',fn = 'amod',
                      MoreArgs = list(mp = c(1,1,0.25 * pi,2),
                                      un = c('mm','m')))
    str(dfm3)
    
    ## seasonal years from 'October' to 'September':
    cl1 <- ringApply(PTclim05,lv = 'year',fn = 'moveYr')
    ##using ringApply to compute multilevel aridity indexes with
    ##function 'wlai':
    tail(cl1,15)
    wl <- ringApply(cl1,lv = 'year',fn = 'wlai')
    str(wl)#only yearly levels with 12 months are evaluated 

    ## A plot of the computed fluctuations of aridity
    d <- groupedData(lmeForm(wl),wl)
    plot(d)        
    
})
