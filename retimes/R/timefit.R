# Stima il picco di una distribuzione creando un istogramma a larghezza fissa:
histestim <- function(x,smoothing=NULL)
{
    if(is.null(smoothing))
        smoothing <- 0
    x <- x[!is.na(x)]
    val <- .C("histestim",
        x = as.double(sort(x)),
        n = as.integer(length(x)),
        smoothing = as.double(smoothing),
        value = as.double(0),
        PACKAGE = "retimes"
    )
    return(val$value)
}

# Stima il picco di una distribuzione valutando le distanze tra ogni x e tutti gli altri
# dati all'interno di un intervallo costruito intorno allo stesso x:
distestim <- function(x,smoothing=NULL)
{
    if(is.null(smoothing))
        smoothing <- 0
    x <- x[!is.na(x)]
    val <- .C("distestim",
        x = as.double(sort(x)),
        n = as.integer(length(x)),
        smoothing = as.double(smoothing),
        value = as.double(0),
        PACKAGE = "retimes"
    )
    return(val$value)
}

# Stima il picco di una distribuzione utilizzando un kernel gaussian:
kernestim <- function(x,smoothing=NULL)
{
    if(is.null(smoothing))
        smoothing <- 0
    x <- x[!is.na(x)]
    val <- .C("kernestim",
        x = as.double(sort(x)),
        n = as.integer(length(x)),
        smoothing = as.double(smoothing),
        value = as.double(0),
        PACKAGE = "retimes"
    )
    return(val$value)
}

# Funzione kernelestim interna:
.kernestim <- function(x,n)
{
    x <- x[!is.na(x)]
    val <- .C("kernestim",
        x = as.double(sort(x)),
        n = as.integer(n),
        h = as.double(0),
        value = as.double(0),
        PACKAGE = "retimes"
    )
    return(val$value)
}

# Metodo dei momenti per l'individuazione dei valori di start per ex-Gaussiana:
mexgauss <- function(x,n=length(x))
{
    k <- start <- c(mu=NaN,sigma=NaN,tau=NaN) # Momenti
    k[1] <- mean(x)
    xdev <- x-k[1]
    k[2] <- sum(xdev^2)/(n-1)
    k[3] <- sum(xdev^3)/(n-1)
    if(k[3] > 0)
        start[3] <- (k[3]/2)^(1/3)
    else
        start[3] <- 0.8*sd(x)
    start[2] <- sqrt(abs(k[2]-start[3]^2))
    start[1] <- k[1]-start[3]
    return(start)
}

# Indica se l'asimmetria è positiva:
skewcheck <- function(x,n)
{
    val <-
        .C("skewcheck",
            x = as.double(x),
            n = as.integer(n),
            check = as.integer(0),
            PACKAGE = "retimes"
        )
    return(val$check)
}

# Calcola l'indice di asimmetria:
skew <- function(x) {
    N <- length(x)
    xdev <- x-mean(x)
    k3 <- sum(xdev^3)/N
    k2 <- sum(xdev^2)/N
    return(k3/(k2^1.5))
}

# Funzione di verosimiglianza:
loglik <- function(par, x)
{
    par[2:3] <- exp(par[2:3])
    D <- (1/par[3])*exp(par[1]/par[3]+(par[2]^2)/(2*par[3]^2)-x/par[3])*
        pnorm(x,par[1]+(1/par[3])*par[2]^2,par[2])
    D <- -sum(log(D))
    if(is.na(D)|is.infinite(D))
        D <- 1e50
    return(D)
}

# Fissati mu e sigma, individua in corrispondenza di quale valore
# di tau si trova il minimo della funzione di verosimiglianza:
loglikMin <- function(x,mu,sigma,tau,n)
{
    L <- rep.int(NA,n)
    for(i in 1:n) {
        D <- dexgauss(x,mu,sigma,tau[i])
        if(sum(is.na(D))==0) {
            if(any(D==0))
                D[which(D==0)] <- 1e-50
            L[i] <- -sum(log(D))
        }
    }
    L[which(L==Inf | L==(-Inf))] <- NA
    tau <- tau[!is.na(L)]
    L <- L[!is.na(L)]
    return(tau[which(L==min(L))][1])
}

# Stima dei parametri distribuzionali per massima verosimiglianza:
timefit <- function(x, iter=0, size=length(x), replace=TRUE, plot=FALSE, start=NULL, ...)
{
    x <- x[!is.na(x)]
    N <- length(x)
    M <- mean(x)
    parNames <- c("mu","sigma","tau")
    if(iter<1) {
        if(is.null(start))
            start <- mexgauss(x,n=N)
        start[2:3] <- log(start[2:3])
        estim <- optim(fn=loglik, par=start, x=x, ...)
        estim$par[2:3] <- exp(estim$par[2:3])
        estim$samples <- estim$bootPar <- matrix()
        estim$bootValue <- NA
        estim$sigmaValid <- NA
    } else {
        start <- c(NA,NA,NA)
        estim <- list(
            par = c(NA,NA,NA),
            value = NA,
            samples = matrix(NA,nrow=size,ncol=iter),
            bootPar = matrix(NA,nrow=iter,ncol=3),
            bootValue = rep.int(NA,iter),
            sigmaValid = NA
        )
        colnames(estim$bootPar) <- parNames
        colnames(estim$samples) <- paste("iter",1:iter,sep="")
        i <- 1
        while(i <= iter) {
            estim$samples[,i] <- sample(x=x,size=size,replace=replace,prob=NULL)
            if(skewcheck(x=estim$samples[,i],n=N)) {
                bootStart <- mexgauss(estim$samples[,i],n=size)
                bootStart[2:3] <- log(bootStart[2:3])
                optimBoot <- optim(fn=loglik, par=bootStart, x=estim$samples[,i], ...)
                estim$bootPar[i,] <- optimBoot$par
                i <- i+1
            }
        }
        # Parameter estimation
        estim$bootPar[,2:3] <- exp(estim$bootPar[,2:3])
        # Mu
        estim$par[1] <- .kernestim(estim$bootPar[,1],n=iter)
        # Sigma
        estim$sigmaValid <- estim$bootPar[,2] >= min(abs(x-M))/sqrt(N-1) & estim$bootPar[,2] <= sd(x)
        if(sum(estim$sigmaValid)<2)
            estim$sigmaValid <- !numeric(iter)
        estim$par[2] <- .kernestim(estim$bootPar[estim$sigmaValid,2],n=sum(estim$sigmaValid))
        # Tau
        estim$par[3] <- loglikMin(x,mu=estim$par[1],sigma=estim$par[2],tau=estim$bootPar[,3],n=iter)
        logParam <- estim$par
        logParam[2:3] <- log(estim$par[2:3])
        estim$value <- loglik(par=logParam,x=x)
        if(plot) {
            par(mfrow=c(2,2))
            for(i in 1:3)
                hist(estim$bootPar[,i],xlab="",main=parNames[i])
        }
    }
    if(plot) {
        D <- density(x)
        lim <- list(x=range(D$x),y=range(D$y))
        lim$y[2] <- lim$y[2]+lim$y[2]/2.5
        plot(D,xlim=lim$x,ylim=lim$y,main="",xlab="",ylab="",lty=1)
        par(new=TRUE)
        curve(dexgauss(x,mu=estim$par[1],sigma=estim$par[2],tau=estim$par[3]),xlim=lim$x,ylim=lim$y,
            main="Distribution",xlab="Observed data",ylab="Density",lty=2)
    }
    names(estim$par) <- parNames
    AIC <- 2*estim$value+2*3
    BIC <- 2*estim$value+3*log(N)
    new("timefit", x=x, samples=estim$samples, par=estim$par, bootPar=estim$bootPar,
        bootValue=estim$bootValue, sigmaValid=estim$sigmaValid, start=start,
        logLik=-estim$value, AIC=AIC, BIC=BIC)
}
