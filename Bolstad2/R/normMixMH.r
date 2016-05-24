normMixMH<-function(theta0, theta1, p,  candidate,
                    steps = 1000, type = 'ind',
                    randomSeed = NULL, startValue = NULL){
    if(steps<100){
        warning("Function should take at least 100 steps")
    }

    if(p<=0 | p>=1)
        stop("Mixture proprotion p must be between 0 and 1")

    mu0<-theta0[1]
    sigma0<-theta0[2]

    mu1<-theta1[1]
    sigma1<-theta1[2]

    mu<-candidate[1]
    sigma<-candidate[2]

    if(any(c(sigma0,sigma0,sigma)<=0))
        stop("All standard deviations must be strictly non-zero and positive")


    if(length(grep('[Ii]',type))>0){
        type <- 'ind'
    }else if(length(grep('[Rr]',type))>0){
        type <- 'rw'
    }else{
        stop("Type must be ind or rw")
    }

    theta<-seq(from = min(mu0-3*sigma0, mu1-3*sigma1),
               to = max(mu0+3*sigma0, mu1+3*sigma1),
               by = 0.001)

    fx<-p*dnorm(theta, mu0, sigma0) + (1-p)*dnorm(theta,mu1,sigma1)
    targetSample<-rep(startValue, steps)

    if(type=='rw'){

        if(!is.null(randomSeed))
            set.seed(randomSeed)

        z<-rnorm(steps, mu, sigma)
        u<-runif(steps)

        if(is.null(startValue))
            startValue <- z[1]

        targetSample[1] <- startValue
        g<-rep(0,steps)
        proposal<-rep(0,steps)
        alpha<-rep(0,steps)

        k1<-p/sigma0*exp(-0.5*((targetSample[1]-mu0)/sigma0)^2)
        k2<-(1-p)/sigma1*exp(-0.5*((targetSample[1]-mu1)/sigma1)^2)
        g[1]<-k1+k2

        i1<-1

        for(n in 2:steps){
            proposal[n]<-targetSample[i1]+z[n]

            k1<-p/sigma0*exp(-0.5*((proposal[n]-mu0)/sigma0)^2)
            k2<-(1-p)/sigma1*exp(-0.5*((proposal[n]-mu1)/sigma1)^2)
            g[n]<-k1+k2

            k3<-g[n]
            k4<-g[i1]

            alpha[n]<-ifelse(k3/k4>1,1,k3/k4)

            ## Metropolis-Hastings Step
            if(u[n]>=alpha[n]){ ## reject
                targetSample[n]<-targetSample[i1]
            }else{ ## accept
                targetSample[n]<-proposal[n]
                i1<-n
            }
        }
    }else{

        if(!is.null(randomSeed))
            set.seed(randomSeed)

        z<-rnorm(steps, mu, sigma)
        u<-runif(steps)

        if(is.null(startValue))
            startValue <- z[1]

        density0<-dnorm(z,mu,sigma)
        density1<-dnorm(z,mu0,sigma0)
        density2<-dnorm(z,mu1,sigma1)

        densityMix<-p*density1+(1-p)*density2

        alpha<-rep(0,steps)
        targetSample[1] <- startValue

        i1<-1
        for(n in 2:steps){
            alpha[n]<-density0[i1]*densityMix[n]/(density0[n]*densityMix[i1])
            alpha[n]<-ifelse(alpha[n]>1,1,alpha[n])

            ## Metropolis-Hastings Step
            if(u[n]>=alpha[n]){
                targetSample[n]<-targetSample[i1]
            }else{
                targetSample[n]<-z[n]
                i1<-n
            }

        }


    }

    oldPar<-par(mfrow=c(1,2),pty="s")

    h<-hist(targetSample, plot = FALSE)
    ymax<-max(c(h$density,fx))*1.05

    hist(targetSample, prob = TRUE, col = "light blue",
         xlim = range(theta), ylim = c(0,ymax),
         main = "Sample from target density",
         xlab = 'x', ylab = 'Density')
    lines(theta, fx)
    box()

    plot(targetSample, type="l", main = "", ylab = "Target Sample")
    par(oldPar)

    invisible(targetSample)
}

