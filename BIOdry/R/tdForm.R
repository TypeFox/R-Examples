tdForm <- structure(function#Time-decline formula
### LME formulas of type time-decline (td) or log-time-decline (ltd).
##details<< (L)td formulas belongs to general equation: log (cs) = log (x) + f(t); where the cummulative radial increments (cs) are explained by the observed radial increments (x) plus a function of time f(t); with f(t) being either the time or a logarithmic transformation the time (ltd).
##references<< Zeide B. 1993. Analysis of Growth Equations. For. Sci., 39: 594-616.                    
(
    rd, ##<<\code{data.frame} object with factor-level columns, or
        ##\code{character} vector with names of factor-level columns
        ##which are decreasingly ordered.
    prim.cov = FALSE, ##<<\code{logical}. Print a primary covariate form:
                      ##'~ cov'. If FALSE then a complete
                      ##formula: 'resp ~ cov | group' is printed.
    on.time = TRUE, ##<< \code{logical}. If TRUE then t =
                    ##'time' (see \code{\link{rtimes}}). If FALSE then
                    ##t = 'year'.
    log.t = FALSE, ##<< \code{logical}. If TRUE then \code{f(t) =
                   ##ln(t)}. Default FALSE produces a td form
    lev.rm = NULL ##<< NULL or \code{character} name of the
                  ##factor-level column to be removed from the
                  ##formula.
) {
    
    rs <- 'log(csx)';lx <- '~ log(x) +' 
    t <- 'time'
    if(!on.time)t <- 'year'
    if(log.t)t <- paste('log(',t,')',sep = '')
    ftt <- paste(rs,lx,t,sep = ' ')
    if(!is.character(rd)){
        f <- colclass(rd,TRUE)[['fac']]
        nf <- rev(f[!f%in%lev.rm])}
    if(is.character(rd))
        nf <- rd[!rd%in%lev.rm]
    if(length(nf) == 0)
        stop('absent levels in formula')
    fc <- paste(nf,collapse = '/')
    fr <- paste(ftt,fc,sep = ' | ')
    if(prim.cov)
        fr <- paste(lx,t,sep = '')
    return(formula (fr))
### \code{formula} with the forms: 'resp ~ cov | group' or '~ cov'.
} , ex=function(){
    ## an ltd formula:
    lev <- c('plot','tree')
    tdeq <- tdForm(lev,log.t = TRUE)
    tdeq
    ## (not run) only primary covariate:
    tdeq1 <- tdForm(lev,prim.cov = TRUE)
    tdeq1
    
    
})
