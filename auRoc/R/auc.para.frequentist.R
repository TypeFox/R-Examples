
auc.para.frequentist.RG <- function(x, y, alpha, method=c("m1", "m2")){

    method <- match.arg(method)

    nx <- length(x)
    ny <- length(y)
    bar.x <- mean(x)
    bar.y <- mean(y)
    s2.x <- var(x)
    s2.y <- var(y)
    hat.M <- (s2.x + s2.y) / (s2.x/nx + s2.y/ny)
    hat.f <- (s2.x + s2.y)^2 / (s2.x^2/(nx-1) + s2.y^2/(ny-1))
    hat.delta <- (bar.x - bar.y) / sqrt(s2.x + s2.y)
    
    ncp <- sqrt(hat.M) * hat.delta
    df <- hat.f

    ci.ncp <- conf.limits.nct(ncp=ncp, df=df, alpha.lower = alpha/2, 
                              alpha.upper = alpha/2)
    ci.delta <- c(hat.delta - qnorm(1-alpha/2)*sqrt(1/hat.M + hat.delta^2/(2*hat.f)),
                  hat.delta + qnorm(1-alpha/2)*sqrt(1/hat.M + hat.delta^2/(2*hat.f)))
    ci <- switch(method,
        m1=c(pnorm(ci.ncp$Lower.Limit/sqrt(hat.M)),
             pnorm(ci.ncp$Upper.Limit/sqrt(hat.M))),
        m2=c(pnorm(ci.delta[1]), pnorm(ci.delta[2])))

    list(PROB=pnorm(hat.delta), C.Interval=as.vector(ci))
}

auc.para.frequentist <- function(x, y, conf.level=0.95, 
                               dist=c("normalDV", "normalEV", "exponential"), 
                               method=c("lrstar", "lr", "wald", "RG1", "RG2")){

    alpha <- 1 - conf.level
    dist <- match.arg(dist)
    method <- match.arg(method)
    dm <- paste(dist, method, sep="")

    estimate <- switch(dm,
                    normalEVwald=Prob(y, x, distr="norm_EV", method="Wald", level=alpha),
                    normalEVlr=Prob(y, x, distr="norm_EV", method="RP", level=alpha),
                    normalEVlrstar=Prob(y, x, distr="norm_EV", method="RPstar", level=alpha),
                    normalDVwald=Prob(y, x, distr="norm_DV", method="Wald", level=alpha),
                    normalDVlr=Prob(y, x, distr="norm_DV", method="RP", level=alpha),
                    normalDVlrstar=Prob(y, x, distr="norm_DV", method="RPstar", level=alpha),
                    normalDVRG1=auc.para.frequentist.RG(x, y, alpha, method="m1"),
                    normalDVRG2=auc.para.frequentist.RG(x, y, alpha, method="m2"),
                    exponentialwald=Prob(y, x, distr="exp", method="Wald", level=alpha),
                    exponentiallr=Prob(y, x, distr="exp", method="RP", level=alpha),
                    exponentiallrstar=Prob(y, x, distr="exp", method="RPstar", level=alpha))

    c(estimate$PROB, estimate$C.Interval)
}
                    

