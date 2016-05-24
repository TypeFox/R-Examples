deviance.lca <-
function(obs, fit, weight, error=c('gaussian', 'poisson'), total=T){
    if (is.character(error)) error <- match.arg(error, c('gaussian', 'poisson'))
    if (error == 'gaussian') res <- weight*(obs-fit)^2
    if (error == 'poisson') {
        ind <- !and(fit > 1e-9, weight) # find cells with fit>0 and weight==1
        fit[ind] <- NA   # set everything else NA so the deviance can be computed:
        res <- 2*weight*(obs*log(obs/fit)-(obs-fit))
    }
    if (total) sum(res, na.rm=T) else res
}
