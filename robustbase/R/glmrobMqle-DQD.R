#### Quasi-Deviance Differences	 --- for Model Selection
#### --------------------------------------------------- -> ./anova-glmrob.R

## MM: These function names are really too long
##     but then, they are hidden in the name space ...

## (Maybe it would be nice to do this as one function with "family" .. )

glmrobMqleDiffQuasiDevB <- function(mu, mu0, y, ni, w.x, phi, tcc)
{
    ##
    f.cnui <- function(u, y, ni, tcc)
    {
	pr <- u/ni
	Vmu <- pr * (1 - pr) ## = binomial()$variance
	residP <- (y-pr)*sqrt(ni/Vmu)

	## First part: nui
	nui <- pmax.int(-tcc, pmin.int(tcc, residP))

	## Second part: Enui
	H <- floor(u - tcc*sqrt(ni*Vmu))
	K <- floor(u + tcc*sqrt(ni*Vmu))
	## Actually, floor is not needed because pbinom() can cope
	## with noninteger values in the argument q!
	## what follows is similar to glmrob.Mqle.EpsiB except a
	## different vectorisation
	h1 <- (if(ni == 1) as.numeric(- (H < 0) + (K >= 1) ) * sqrt(Vmu)
	       else
	       (pbinom(K-1,1,pr) - pbinom(H-1,ni-1,pr)
		- pbinom(K,ni,pr) + pbinom(H,ni,pr)) * pr * sqrt(ni/Vmu))
	## pmax was needed to get numeric returns from pbinom

	Enui <- (tcc*(1 - pbinom(K,ni,pr) - pbinom(H,ni,pr)) + h1)

	return((nui - Enui) / sqrt(ni*Vmu))
    } ## f.cnui()

    nobs <- length(mu)
    stopifnot(nobs > 0)
    QMi <- numeric(nobs)
    ## Numerical integrations
    for(i in 1:nobs)
	QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], tcc = tcc,
			    subdivisions = 200,
			    lower = mu[i]*ni[i], upper = mu0[i]*ni[i])$value
    ## robust quasi-deviance
    ## -2*(sum(QMi1)-sum(QMi2))	   ## Andreas' interpretation of (4) and (5)
    ## -2*(sum(QMi1)-sum(QMi2)/nobs)  ## Eva's interpretation of (4) and (5)
    ## According to Andreas' interpretation
    -2*sum(QMi*w.x)
} ## glmrobMqleDiffQuasiDevB

glmrobMqleDiffQuasiDevPois <- function(mu, mu0, y, ni, w.x, phi, tcc)
{
    ##
    f.cnui <- function(u, y, ni, tcc)
    {
	Vmu <- u ## = poisson()$variance
	residP <- (y-u)/sqrt(Vmu)

	## First part: nui
	nui <- pmax.int(-tcc, pmin.int(tcc, residP))

	## Second part: Enui
	H <- floor(u - tcc*sqrt(Vmu))
	K <- floor(u + tcc*sqrt(Vmu))
	## what follows is similar to Epsipois except a
	## different vectorisation
	h1 <- u/sqrt(Vmu)*(dpois(H,u)- dpois(K,u))
	Enui <- tcc*(1 - ppois(K,u) - ppois(H,u)) + h1

	return((nui - Enui) / sqrt(Vmu))
    }

    nobs <- length(mu)
    stopifnot(nobs > 0)
    QMi <- numeric(nobs)
    ## Numerical integrations
    for(i in 1:nobs)
	QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], tcc = tcc,
			    lower = mu[i], upper = mu0[i])$value

    ## robust quasi-deviance
    ## -2*(sum(QMi1)-sum(QMi2))	  ## Andreas' interpretation of (4) and (5)
    ## -2*(sum(QMi1)-sum(QMi2)/nobs) ## Eva's interpretation of (4) and (5)
    ## According to Andreas' interpretation
    -2*sum(QMi*w.x)
}## glmrobMqleDiffQuasiDevPois

glmrobMqleDiffQuasiDevGamma <- function(mu, mu0, y, ni, w.x, phi, tcc,
                                        variant = c("V1", "Eva1", "Andreas1"))
{
    ## Notation similar to the discrete case (Cantoni & Ronchetti, 2001)
    f.cnui <- function(u, y, ni, phi, tcc)
    {
        s.ph <- sqrt(phi)
	## First part: nui
	sV <- s.ph * u ## = sqrt(dispersion * Gamma()$variance)
	residP <- (y-u)/sV
	nui <- pmax.int(-tcc, pmin.int(tcc, residP))

	## Second part: Enui
        ## what follows is similar to glmrob.Mqle.Epsipois except a
	## different vectorisation
        nu <- 1/phi      ## form parameter nu
        snu <- 1/s.ph ## sqrt (nu)

        pPtmc <- pgamma(snu - tcc, shape=nu, rate=snu)
        pPtpc <- pgamma(snu + tcc, shape=nu, rate=snu)

        Enui <- tcc*(1-pPtpc-pPtmc) + Gmn(-tcc,nu) - Gmn( tcc,nu)

        ( nui/sV - Enui/u*s.ph )
    }
    f.cnui1 <- function(u, y, ni, phi, tcc)
    {
	## First part: nui
	sV <- sqrt(phi) * u ## = sqrt(dispersion * Gamma()$variance)
	residP <- (y-u)/sV
	nui <- pmax.int(-tcc, pmin.int(tcc, residP))
        (nui  / sV)
    }
    f.cnui2 <- function(u, y, ni, phi, tcc)
    {
	## First part: nui
        s.ph <- sqrt(phi)
	sV <- s.ph * u ## = sqrt(dispersion * Gamma()$variance)
        snu <- 1/s.ph ## sqrt (nu)

	## Second part: Enui
        ## what follows is similar to EpsiGamma except a
	## different vectorisation
        nu <- 1/phi ## form parameter nu

        pPtmc <- pgamma(snu - tcc, shape=nu, rate=snu)
        pPtpc <- pgamma(snu + tcc, shape=nu, rate=snu)

        Enui <- tcc*(1-pPtpc-pPtmc) + Gmn(-tcc,nu) - Gmn( tcc,nu)
	return(Enui/u * s.ph)
    }

    nobs <- length(mu)
    stopifnot(nobs > 0)
    variant <- match.arg(variant)
    ## robust quasi-deviance
    if(variant == "V1") {
        QMi <- numeric(nobs)
        ## Numerical integrations
        for(i in 1:nobs)
            QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
                                lower = mu[i], upper = mu0[i])$value
        -2*sum(QMi*w.x)

    } else { ## "Eva1" or "Andreas1";  Using two terms
        QMi1 <- QMi2 <- numeric(nobs)
        for(i in 1:nobs)
            QMi1[i] <- integrate(f.cnui1, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
                                 lower = mu[i], upper = mu0[i])$value
        for(i in 1:nobs)
            QM2i[i] <- integrate(f.cnui2, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
                                 lower = mu[i], upper = mu0[i])$value

        if(variant == "Eva1") {  ## Eva Cantoni's interpretation of (4) and (5)
            -2*(sum(QMi1)-sum(QMi2)/nobs)
        } else if (variant == "Andreas1") { ## Andreas' interpretation of (4) and (5)
            -2*(sum(QMi1)-sum(QMi2))
        } else stop("invalid 'variant': ", variant)
    }
}

