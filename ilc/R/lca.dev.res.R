lca.dev.res <-
function(lca.obj, pop, clip=0){
    obs <- pop*lca.obj[[4]]
    fit <- lca.obj$fit$y
    if (any(fit<0)) fit <- exp(fit)
    if (any(dim(pop) != dim(fit))) stop('Wrong population data!')
    fit <- pop*fit
    mw <- bool(pop > 1e-9, na=F)
    if (clip) { 
        year <- lca.obj$year; n <- length(year)
        age <- lca.obj$age; k <- length(age)
        ft <- gl(n, k)
        stx <- year[1]-age[k]   # year of birth of first cohort  
        etx <- year[n]-age[1]   # year of birth of last cohort 
        mtx <- array(year[ft]-age, c(k,n)) #  matrix(year[ft]- age[fx], ncol=n, nrow=k)
        mw[bool(mtx < stx+clip)] <- F   
        mw[bool(mtx > etx-clip)] <- F
        warning(paste('The cohorts outside',  sqb(paste(c(stx+clip, etx-clip), 
                                                        collapse=', ')), 'were zero weighted (clipped).'))
    }
    dev <- deviance.lca(obs, fit, mw, error='poisson', total=FALSE)    
    lca.obj$residuals$y <- sign(obs-fit)*sqrt(dev/lca.obj$mdev[1])
    lca.obj
}
