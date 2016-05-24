singlePerM <- function(tt,yy,s,n,ii, pp,W,rho, regression, design, steps, var1, tol, genoudcontrol) {
    # Used by RobPer 
    # Calculates one periodogram bar for given trial period pp
    # Arguments not checked as it is assumed that the function RobPer gives the right arguments

    # Initial Full model:
    # Build design matrix
    XX <- Xgen(tt=tt,n=n,s=s, pp=pp, design=design, steps=steps)
    if(is.na(XX[1])) {
        return(NA)
    } else {
    	
        # Initial LTS, small nsamp in order to be faster
        if(all(dim(XX)%%2==c(1,0))) alpha <- 0.5-1/(n-dim(XX)[2]-1) else alpha <- 0.5
        tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) #ltsReg seldomly crashes... 
        if(inherits(tempFULL, "try-error"))  tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) #give it a second chance
        if(inherits(tempFULL, "try-error"))  tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) # and a third time
        if(inherits(tempFULL, "try-error"))  tempFULL <- suppressWarnings(rq(yy~0+XX, tau=.5, method="br")) ## if still problems, use L1
    	
        rr <- as.vector(tempFULL$residuals)
    	
        # Estimate of the standard deviation
        if(var1) sigma <- 1 else sigma <- median(abs(rr[rr!=0]))/0.675
    	
        # Initial L1 regression of the intercept model
        tempINT <- suppressWarnings(rq(yy~0+ii, tau=.5, method="br"))
        ee <- as.vector(tempINT$residuals)
    	
        # IRWLS for intercept model
    	
        if (regression=="bisquare") {
            gamma <- lmrob..M..fit(x=ii, y=yy, beta.initial=tempINT$coeff, scale=sigma,
                control=lmrob.control(tuning.psi=4.685061, psi="bisquare"))$coeff
        }
        else {
            gamma <- IRWLS(yy, matrix_=ii, W=W, residuals_=ee, scale_=sigma, tol=tol)
        }
    	
        # Full model continueing
        # Use genetic optimization in case of biweight regression in order to get a good starting point
        if(regression=="bisquare") {
            beta_ <- tempFULL$coeff
            beta_gamma <- numeric(length(beta_))
            if(design %in% c("step", "stepB","splines")) beta_gamma<- beta_gamma+gamma
            if(design %in% c("sine", "fourier(2)", "fourier(3)")) beta_gamma[1]<- gamma
            optf <- function(beta) sum(rho((yy-XX%*%beta)/sigma))
            beta_ <- suppressWarnings(genoud(optf, nvars=length(beta_),starting.values=rbind(beta_, beta_gamma),
                pop.size=genoudcontrol$pop.size, max=FALSE, max.generations=genoudcontrol$max.generations,
                wait.generations=genoudcontrol$wait.generations, solution.tolerance=tol, print.level=0)$par)
            rr<- as.vector(yy-XX%*%beta_)
        }
    	
        # IRWLS for full model
        if (regression=="bisquare") {
            beta <- lmrob..M..fit(x=XX, y=yy, beta.initial=beta_, scale=sigma, control=lmrob.control(tuning.psi=4.685061, psi="bisquare"))$coeff
        }
        else {
            beta  <- IRWLS(yy, matrix_=XX, W=W, residuals_=rr, scale_=sigma, tol=tol)
        }
    	
        bar<- 1- sum(rho((yy-XX%*%beta)/sigma))/ sum(rho((yy-ii*gamma)/sigma))
    	
        return(bar)
    }
}
