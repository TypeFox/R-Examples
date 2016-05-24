lca.set <- function(dat, year=dat$year, age=dat$age, series=1, 
                    max.age=100, interpolate = F){
    
    # check data type:
    if (class(dat) != "demogdata" || dat$type != "mortality")
        stop("Not mortality data")
    
    error <- get('error', pos=parent.frame())
    
    cat('Original sample: ')
    print(dat)
    
    # check data ranges:
    ind <- year%in%dat$year
    if (!all(ind)) warning(paste('The following years are not available:',
                                 paste(year[!ind], collapse=' '), sep='\n'))
    ind <- age%in%dat$age
    if (!all(ind)) warning(paste('The following ages are not available:',
                                 paste(age[!ind], collapse=' '), sep='\n'))
    gr <- F
    if (max.age<max(age)) {
        age <- min(age):max.age; gr <- T
        warning(paste(' => data above age', max.age, 'are grouped.'))
    }
    
    # Restrict data set to the given age/period.
    # Note: the (demography) function below extracts all available series in dat
    # wich is not really necessary since only one series is fitted!
    # However, I've kept it for now to make things simpler.
    if (!all(dat$year%in%year)) dat <- extract.years(dat, year)
    year <- dat$year
    if (!all(dat$age%in%age)) dat <- extract.ages(dat, age, combine.upper=gr)
    age <- dat$age
    
    # first extract and check the force of mortality:
    mx <- dat$rate[[series]]
    ind <- bool(mx < 1e-09, na=T) # get 0 (or negative) and NA value index
    if (sum(ind)) { 
        if (interpolate) {
            # re-estimate all 0 and NA mu-s:
            dat <- fill.demogdata(dat, series=series, method='interp') # mspline
            mx <- dat$rate[[series]]
            ind <- bool(mx < 1e-09, na=T) # reset index for final check
            if (sum(ind)) stop('It seems that interpolation did not help!\n')
        } else {
            # re-set all 0 mu-s to NA (due to log(0)=Inf )
            mx[ind] <- NA
            warning(paste("There are", sum(ind), "cells with 0/NA mu, which are", 
                          "ignored in the current analysis.\n  Try reducing the maximum age", 
                          "or setting interpolate=TRUE."))
        }
    }    
    mz <- log(mx)
    
    # Once the mu-s are checked, extract the exposure matrix and the weights:
    me <- dat$pop[[series]]
    mw <- if (error == 'poisson') bool(me>1e-9, na=F) else array(1, dim=dim(mx)) 
    # for poisson errors: weights = 0/1 <=> exp == 0
    # for gaussian errors: weights = 1 (except for NA mu's)
    if (sum(!mw)){
        warning(paste("There are", sum(!mw), "cells with 0/NA exposures,", 
                      "which are ignored in the current analysis.\n  Try reducing the", 
                      "maximum age or choosing a different age range.\n  Alternatively, fit LC", 
                      "model with error=", mark('gaussian', F), "."))
    }
    
    cat('Applied sample: ')
    print(dat)
    
    # Finally, get the number of deaths
    # insp.dd() will make sure to get 0 deaths for 0 exposures (i.e. before this
    # was re-estimated in the case of interpolate)
    md <- insp.dd(dat, what='d', series=series)
    n <- length(year); k <- length(age); rtx <- n+k-1  # cohort range
    assign('mx', mx, pos=sys.frame(-1))
    assign('mz', mz, pos=sys.frame(-1))
    assign('me', me, pos=sys.frame(-1))
    assign('md', md, pos=sys.frame(-1))
    assign('mw', mw, pos=sys.frame(-1))
    assign('age', age, pos=sys.frame(-1))
    assign('year', year, pos=sys.frame(-1))
    assign('n', n, pos=sys.frame(-1)) 
    assign('k', k, pos=sys.frame(-1))
    assign('rtx', rtx, pos=sys.frame(-1))
    return(list(n = n, k = k, rtx = rtx, mx = mx, md = md, mz = mz, me = me))
}
