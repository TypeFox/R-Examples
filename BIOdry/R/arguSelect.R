arguSelect <- structure(function#Argument selection
###Getting arguments in function(s).                        
##details<< This function is implemented by \code{\link{ringApply}} to model multilevel data in \code{rd}. But, the function may also be used for other purposes.  
(
    rd = NULL, ##<<\code{NULL} or \code{data.frame} object with
               ##factor-level columns (see details).
    fun = c('mapply','ringApply'), ##<<\code{character} vector with
    ## name(s) of formulated functions.
    ... ##<< Further arguments not necessarily required by formulated
        ##function(s). Arguments in \code{MoreArgs} lists are also
        ##processed.
) {
    mx <- list(...)
    
    if('y'%in%names(mx)){
        nl <- names(splitFrame(rd))
        refs <- levexp(mx[['y']],nl)
        mx[['y']] <- refs[nl]}
    
    mar <- 'MoreArgs'
    fn <- mx[['fn']]
    fun <- c(fun,fn)
    fca <- lapply(fun,
                  function(x)names(formals(x)))
    nfr <- unlist(fca)
    mx. <- mx[!names(mx)%in%mar]
    sel <- mx.[names(mx.)%in%nfr]
    s <- names(mx[[mar]])%in%nfr
    if(any(s))
        sel[[mar]] <- mx[[mar]][s]
    if(is.data.frame(rd))
        sel[['rd']] <- rd
    return(sel)
### \code{list} of arguments.
} , ex=function() {
    
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    ## Radial increments measured on 2003:
    data(Pradii03,envir = environment())    

    ## getting arguments in some functions:
    ar1 <- arguSelect(fun = c('amod'),
                      only.dup = TRUE,mp = c(0.5,1),z = 2003)
    str(ar1)
    
    ar2 <- arguSelect(fn = 'amod',
                      only.dup = TRUE,mp = c(0.5,1),z = 2003)
    str(ar2)
    ar3 <- arguSelect(rd = Prings05,fn = 'amod',
                      only.dup = TRUE,mp = c(0.5,1),z = 2003)
    str(ar3)
    
    ar4 <- arguSelect(rd = Prings05,
                      fun = 'scacum',y = Pradii03,
                      MoreArgs = list(only.dup = TRUE,
                                      mp = c(0.5,1),z = 2003))
    str(ar4)
    
    ar5 <- arguSelect(rd = Prings05,
                      fun = 'scacum',y = Pradii03,z = rep(2003:2011),
                      MoreArgs = list(only.dup = TRUE,
                                      mp = c(0.5,1)))
    str(ar5)    
    
})
