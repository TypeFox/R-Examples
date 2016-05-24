rhdata <-
function(dat, covar, xbreaks=60:96, xlabels=60:95,
                   ybreaks=mdy.date(1,1,1999:2008), ylabels=1999:2007,  
                   name=NULL, label=NULL){
    
    # make sure all covariates are factors
    if (mode(covar) == 'character'){
        tmp <- pmatch(covar, names(dat))
        ind.na <- is.na(tmp)
        if (any(ind.na)) 
            stop(paste('Variates', paste(covar[ind.na], collapse=', '), 
                       'could not be found!'))
        covar <- tmp
    } 
    # make sure that all covariates are factors
    # Note: it's not efficient to transform the variates into factors here
    # as some of them could be continuous or characters, etc.
    for(i in covar)
        if (!is.factor(dat[[i]])) stop('All covariates must be factors!\n')
    
    # browser()
    
    # calculate central exposure:
    tmp <- rep(0, nrow(dat))
    ysc <- 365.25 # year scale
    cut.age <- tcut(tmp, floor(xbreaks * ysc), labels=xlabels)
    cut.yr <- tcut(dat$dob, ybreaks, labels=ylabels)
    fu <- with(dat, dev+1-dob)    # follow up time 
    # Note: the +1 ensures that the date of death (event) is included as a whole day
    expr <- paste('fu~cut.age+cut.yr',paste(names(dat)[covar], collapse='+'), sep='+')
    e.xy <- with(dat, pyears(as.formula(expr), scale=ysc)$pyears)
    
    # calculate age-year-specific number of deaths: only for T events
    # for safety remove NAs, though there shouldn't be any in the data sets
    ind.na <- bool(dat$event, na=F)
    # detach(); attach(dat[ind,])
    fu[!ind.na] <- 0   # make the follow up times 0 for censored cases
    cut.age <- tcut(fu, floor(xbreaks * ysc), labels=xlabels)
    cut.yr <- tcut(dat$dev, ybreaks, labels=ylabels) 
    tmp[ind.na] <- 1 # create a fictive follow up time of 1 day  
    fu <- tmp # reset fu so we can keep the same expression
    n.xy <- with(dat, pyears(as.formula(expr), scale=ysc)$n)
    mu  <-  n.xy/e.xy
    mu <- ifelse(is.infinite(mu), NA, mu)
    # browser()
    
    # return data object 
    ret <- list(year=ylabels, age=xlabels, covariates = dimnames(mu)[-(1:2)], 
                deaths=n.xy, pop=e.xy, mu=mu, label = label)
    if (!is.null(name)) ret$name <- name
    names(ret$covariates) <- names(dat)[covar]
    return(structure(ret, class='rhdata'))
}
