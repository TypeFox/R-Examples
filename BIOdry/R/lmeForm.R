lmeForm <- structure(function#LME formula
### LME formula with grouping levels being defined by factor-level
### columns in a data frame.
##references<< Pinheiro J. C., D. M. Bates. 2000. Mixed-effects models in S and S-PLUS. Springer, New York.                     
(
    rd, ##<< \code{data.frame} object with fator-level columns.
    prim.cov = FALSE, ##<<\code{Logical}: should the formula be
                      ##printed in the form of primary covariate: '~
                      ##cov'?  If FALSE then a complete form: 'resp
                      ##~ covar | group' is formulated.
    resp = NULL, ##<<\code{NULL} or \code{character}. Column name of
                 ##the response. If NULL name of the first numeric
                 ##column of the data frame is used.
    covar = NULL, ##<<\code{NULL} or \code{character}. Column name(s)
                  ##of the covariate(s). If \code{NULL} name of the
                  ##first temporal column in the data frame is used.
    lev.rm = NULL ##<< \code{NULL}, \code{character} or \code{numeric}
                  ##vector of levels in the data frame to be removed
                  ##from the groups.
    
) {
    if(is.null(resp))
        resp <- colclass(rd)[['num']][1]
    if(is.null(covar))
        covar <- colclass(rd)[['tmp']][1]        
    covar. <- paste('~',covar,sep = ' ')
    covar <- paste(resp,'~',covar,sep = ' ')
    f <- colclass(rd)[['fac']]
    if(is.numeric(lev.rm))
        lev.rm <- f[lev.rm]
    nf <- rev(f[!f%in%lev.rm])
    if(length(nf) == 0)
        stop('absent levels in formula')
    fc <- paste(nf,collapse = '/')
    fr <- paste(covar,fc,sep = ' | ')
    fr <- formula(fr,showEnv = FALSE)
    if(prim.cov)fr <- covar.

    return(fr)
### \code{formula} with any of the forms: \code{resp ~ cov | group} or
### \code{~ cov}.
} , ex=function(){
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,,envir = environment())
    
     ## groupedData object with  lmeForm 
    gdata <- groupedData(lmeForm(Prings05,lev.rm = 1),data = Prings05)
    plot(gdata,groups = ~ sample)

    ## LME formula in form of covariate:
    form1 <- lmeForm(Prings05,prim.cov = FALSE)
    print(form1)

})
