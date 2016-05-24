
Lognorm.diff<-function(x,y, conf.level=0.95, alternative="two.sided", sim=10000,...)
{
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
    args <- list(...)
    nx <- length(x)
    ny <- length(y)

    if (all(x > 0) & all(y > 0)) {
        lx <- log(x)
        ly <- log(y)
        active <- TRUE
    }
    else {
        if (any(x < 0) | any(y < 0)) {
            estimate <- NA
            conf.int <- c(NA, NA)
            active = FALSE
            warning("negative values occured")
        }
        else {
            lx <- log(x + 0.1)
            ly <- log(y + 0.1)
            active <- TRUE
            warning("0.1 added to x and y, because 0 accured")
        }
    }

if (active)
    {
    mlx <- mean(lx)
    varlx <- var(lx)
    mly <- mean(ly)
    varly <- var(ly)
    mx <- exp(mlx + 0.5 * varlx)
    my <- exp(mly + 0.5 * varly)
    estimate <- mx - my

    Zx<-rnorm(n=sim, mean=0, sd=1)
    Chix<-rchisq(n=sim, df=nx-1)
    Tx <- mlx - (Zx*sqrt(varlx))/((sqrt(Chix)/sqrt(nx-1))*sqrt(nx)) + (varlx)/(2*Chix/(nx-1))

    Zy<-rnorm(n=sim, mean=0, sd=1)
    Chiy<-rchisq(n=sim, df=ny-1)
    Ty <- mly - (Zy*sqrt(varly))/((sqrt(Chiy)/sqrt(ny-1))*sqrt(ny)) + (varly)/(2*Chiy/(ny-1))

    TD <- exp(Tx)-exp(Ty)

    switch(alternative,
    two.sided={conf.int<-quantile(x=TD, probs=c((1-conf.level)/2, 1-(1-conf.level)/2))},
    less={conf.int<-c( -Inf, quantile(x=TD, probs=conf.level) )},
    greater={conf.int<-c( quantile(x=TD, probs=1-conf.level), Inf )}
    )
    }

METHOD<-"Difference of means assuming lognormal distribution"

attr(conf.int, which="methodname")<-METHOD

    return(list(conf.int = conf.int, estimate = estimate))
}