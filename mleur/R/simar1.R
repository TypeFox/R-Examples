simar1 <-
function (phi = 0.5, n = 100, InnovationVariance = 1, 
            noiseDist=c("normal", "t", "stable", "GARCH11"), 
            df=5, ALPHA=1.5, BETA=0, GAMMA=1, DELTA=0, alpha=0.2, beta=0.7) 
{
    p <- length(phi)
    if (length(p) != 1) stop("error: invalid phi")
    noiseD <- noiseDist[match(noiseDist, noiseDist)[1]]
    sdA <- sqrt(InnovationVariance)
#generate white noise
    spec <- garchSpec(model = list(alpha = alpha, beta = beta))
    a  <- switch(noiseD,
            normal = rnorm(n,mean = 0, sd = sqrt(InnovationVariance)),
            t = rt(n=n, df=df),
            stable = rstable(n, ALPHA, BETA,GAMMA,DELTA),
            GARCH11 =as.vector(garchSim(spec, n = n))
            )
#random walk case
    if (abs(phi) >= 1) 
        return(cumsum(a))
#get z[1], the initial value
#This is for the non-normal case, normal case below.
#First initial values to innovations.
    Q <- 2^(ceiling(log(n, base=2))+1)-n
    psi <- phi^(0:(Q-1))
    a1 <- switch(noiseD,
            normal = NA,
            t = rt(n=Q+1, df=df),
            stable = rstable(Q+1, ALPHA, BETA),
            GARCH11 =as.vector(garchSim(spec, n = Q+1))
            )
#now z[1]
    if (noiseD=="normal") z1<-rnorm(1)*sdA/sqrt(1-phi^2) else 
        z1 <- (rev(convolve(a1, rev(psi), type = "o"))[-(1:Q)])[1]
    z <- numeric(n)
    z[1] <- z1
#remaining z's
    for (i in 2:n) z[i] = a[i] + phi * z[i - 1]
    z
}

