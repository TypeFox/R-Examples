"km.ci" <-
function(survi,conf.level=0.95, tl=NA, tu=NA, method="rothman")
{
    # This function can compute the most desirable confidence bands.
    # The method "log" is implemented as "log" in R survfit.
    # The method "loglog" is implemented as "log-log" in R survfit.
    # The method "linear" is called "plain" in R survfit.

    if(conf.level < 0 || conf.level > 1)
      stop("confidence level must be between 0 and 1")
    if (data.class(survi)!="survfit")
            stop("Survi must be a survival object")
    method <- match.arg(method,c( "peto", "linear", "log" ,"loglog", "rothman","grunkemeier",
                "epband", "logep", "hall-wellner","loghall"))

    if(method=="grunkemeier")
    {
        result <- grunk.all.fun(survi,1-conf.level)
        result$conf.type <- "Grunkemeier"
    }

    if(method=="linear")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$linear$lower
        result$upper <- cf$linear$upper
        result$conf.type <- "Linear"
    }

    if(method=="rothman")
    {
        result <- rothman.fun(survi,conf.level)$surv.object
        result$conf.type <- "Rothman"
    }
    if(method=="peto")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$peto$lower
        result$upper <- cf$peto$upper
        result$conf.type <- "Peto"
    }
    if(method=="log")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$greenwood$lower
        result$upper <- cf$greenwood$upper
        result$conf.type <- "Log"
    }
    if(method=="loglog")
    {
        result <- survi
        cf <- comp.npci(survi,conf.level)
        result$lower <- cf$log$lower
        result$upper <- cf$log$upper
        result$conf.type <- "Log-Log"
    }
    if(method=="hall-wellner")
    {
        result <- hall.wellner.fun(survi, tl=tl, tu=tu, conf.lev=conf.level)
        result$conf.type <- "Hall-Wellner"
    }
    if(method=="loghall")
    {
        result <- hall.wellner.fun(survi,tl=tl, tu=tu, method="log", conf.lev=conf.level)
        result$conf.type <- "Log(Hall-Wellner)"
    }

    if(method=="epband")
    {
        result <- epband.fun(survi, tl=tl, tu=tu, conf.lev=conf.level)
        result$conf.type <- "Equal Precision"
    }
    if(method=="logep")
    {
        result <- epband.fun(survi, tl=tl, tu=tu, method="log",conf.lev=conf.level)
        result$conf.type <- "Log(Equal Precision)"
    }
    return(result)
}

