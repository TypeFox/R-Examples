`icfitControl`<-function (epsilon = 1e-06, maxit = 10000, initfitOpts=NULL, 
      conf.level=.95, B=200, confMethod="modboot", seed=19439101, timeEpsilon=1e-06, timeMessage=TRUE) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(conf.level) & (conf.level>=1 | conf.level<=0))
        stop("conf.level must be between 0 and 1")
    if (!is.numeric(B)  || B<=10)
        stop("B must be at least 11")

    list(epsilon = epsilon, maxit = maxit, initfitOpts= initfitOpts, conf.level=conf.level,
        B=B, confMethod=confMethod, seed=seed, timeEpsilon=timeEpsilon, timeMessage=timeMessage)
}
