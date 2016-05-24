ruinprob.test <- function(x, prob.null, type = c("bootstrap", "normal"), nboot, bootmethod = c("nonp", "exp", "lnorm"), ...){
    ruinargs <- list(...)
    type <- match.arg(type)

    if(length(ruinargs$se.multiplicator)){
        semult   <- ruinargs$se.multiplicator
        ruinargs <- ruinargs[-which(names(ruinargs) == "se.multiplicator")]
    } else {
        semult <- 4.0
    }

    jackint  <- semult * ifelse(exists("ruinargs$interval"),
                                ruinargs$param.list$interval,
                                ruinargs$interval)
    dataname <- deparse(substitute(x))

    #x        <- rpdataconv(x)
    rp       <- ruinprob(x = x, param.list = ruinargs)
    rp.se    <- rpjack(x = x, param.list = ruinargs, interval = jackint)

    teststat <- rpteststat(rp, prob.null, rp.se)

    if(type == "bootstrap"){
        x.boot        <- rpdataboot(x = x, b = nboot, method = ifelse(missing("bootmethod"), "nonp", bootmethod))
        rp.boot       <- ruinprob(x = x.boot, param.list = ruinargs)
        rp.boot.se    <- rpjack(x = x.boot, param.list = ruinargs, interval = jackint)
        teststat.boot <- rpteststat(rp.boot, rp, rp.boot.se)
    }

    structure(
        list(statistic   = structure(teststat, names = "T"),
             parameter   = unlist(ruinargs[vapply(ruinargs, is.numeric, logical(1L))]),
             p.value     = switch(type,
                  normal    = rppvalue(teststat, "normal"),
                  bootstrap = rppvalue(teststat, "bootstrap", teststat.boot)),
             estimate    = structure(rp, names = "probability of ruin"),
             null.value  = structure(prob.null, names = "probability of ruin"),
             alternative = "less",
             method      = paste("A test for the probability of ruin using", type, "approximation."),
             data.name   = dataname
        ),
        class = "htest"
    )
}
