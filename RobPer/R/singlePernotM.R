singlePernotM <- function(tt,yy,s,n,ii, pp, regression, design, steps, zeta, LTSopt, tol,  genoudcontrol,  gamma, K_gamma, taucontrol, Scontrol) {
    # Used by RobPer
    # Calculates one periodogram bar for given trial period pp
    # when performing other regression than M regression
    # Arguments not checked as it is assumed that the function RobPer gives the right arguments
    
    # Build design matrix
    XX <- Xgen(tt=tt,n=n,s=s, pp=pp, design=design, steps=steps)
    if(is.na(XX[1])) {
        return(NA)
    }
    if(dim(XX)[2]==1) {
        return(0)  # In case of only one step to fit
    }
    if(!is.na(XX[1])&(dim(XX)[2]!=1)) {
    	
        if(regression=="LTS"){
            # Intercept model
            if(is.na(gamma)) { # in case of step model
                alpha <- (n+dim(XX)[2]-3)/(2*(n-2)) # use trimming h of full model
                if(length(unique(ii))==1) {
                    tempINT  <- ltsReg(yy~1  , use.correction=FALSE, alpha=alpha)
                } else {
                    tempINT  <- ltsReg(yy~0+ii,use.correction=FALSE, alpha=alpha)
                }
                gamma <- tempINT$coeff
                ee <- tempINT$residuals
                K_gamma <- zeta(ee,  q=dim(XX)[2])
            }
    		
            # Full model
            LTSproblem <- FALSE
            if(all(dim(XX)%%2==c(1,0))) alpha <- 0.5-1/(n-dim(XX)[2]-1) else alpha <- 0.5
            if(all(XX[,1]==1)) {
                tempFULL<-try(tempFULL<-ltsReg(yy~1+XX[,-1], mcd=FALSE, alpha=alpha), silent=TRUE) #ltsReg seldomly crashes... give it a second chance
                if(inherits(tempFULL, "try-error")) {
                    tempFULL<-try(ltsReg(yy~1+XX[,-1], mcd=FALSE, alpha=alpha), silent=TRUE)
                } #... and a third time
                if(inherits(tempFULL, "try-error")) {
                    tempFULL<-try(ltsReg(yy~1+XX[,-1], mcd=FALSE, alpha=alpha), silent=TRUE)
                }
                if(inherits(tempFULL, "try-error")) {
                    LTSproblem<-TRUE
                } # if still problems, skip this period.
            } else {
                tempFULL<-try(ltsReg(yy~0+XX, mcd=FALSE, alpha=alpha), silent=TRUE) #ltsReg seldomly crashes... give it a second chance
                if(inherits(tempFULL, "try-error")) {
                    tempFULL<-try(ltsReg(yy~0+XX, mcd=FALSE, alpha=alpha), silent=TRUE)
                } #... and a third time
                if(inherits(tempFULL, "try-error")) {
                    tempFULL<-try(ltsReg(yy~0+XX, mcd=FALSE, alpha=alpha), silent=TRUE)
                }
                if(inherits(tempFULL, "try-error")){  LTSproblem<-TRUE} # if still problems, skip this period.
            }
            if(LTSproblem) {
                if(LTSopt) {
                    tempFULL<-suppressWarnings(rq(yy~0+XX, tau=.5, method="br"))
                } else {
                    return(NA)
                    warning(paste("LTS-fit did not work for trial period", pp))
                }
            }
    		
            # Use genetic optimization
            if(LTSopt) {
                beta_ <- tempFULL$coeff
                beta_gamma <- numeric(length(beta_))
                if(design %in% c("step", "stepB","splines")) beta_gamma<- beta_gamma+gamma
                if(design %in% c("sine", "fourier(2)", "fourier(3)")) beta_gamma[1]<- gamma
                optf <- function(beta) zeta(yy-XX%*%beta, dim(XX)[2])
                beta_ <-suppressWarnings(genoud(optf, BFGS=FALSE, nvars=length(beta_),starting.values=rbind(beta_,beta_gamma ),
                    pop.size=genoudcontrol$pop.size, max=FALSE, max.generations=genoudcontrol$max.generations,
                    wait.generations=genoudcontrol$wait.generations, solution.tolerance=tol, print.level=0)$par)
                rr<- as.vector(yy-XX%*%beta_)
            } else rr<- tempFULL$residuals
            K_beta<- zeta(rr,  q=dim(XX)[2])
        }
    	
        if(regression=="L2") {
            # Full model
            rr <- lm(yy~0+XX)$residuals
            K_beta <- zeta(rr)
        }
    	
        if(regression=="L1") {
            rr <- suppressWarnings(rq(yy~0+XX, tau=.5, method="br")$residuals)
            K_beta<- zeta(rr)
        }
    	
        if(regression=="S") {
            beta_gamma <- numeric(dim(XX)[2])
            if(design %in% c("step", "stepB","splines")) beta_gamma<- beta_gamma+gamma
            if(design %in% c("sine", "fourier(2)", "fourier(3)")) beta_gamma[1]<- gamma
            K_beta<-FastS(XX,yy,Scontrol=Scontrol, beta_gamma=beta_gamma)$scale
        }
    	
        if(regression=="tau"){
            beta_gamma <- numeric(dim(XX)[2])
            if(design %in% c("step", "stepB","splines")) beta_gamma<- beta_gamma+gamma
            if(design %in% c("sine", "fourier(2)", "fourier(3)")) beta_gamma[1]<- gamma
            K_beta <- FastTau(XX, yy, beta_gamma=beta_gamma, taucontrol=taucontrol)$scale
            if(is.na(K_beta))  warning(paste("trial period", pp, "could not be fitted, too many degenerate subsamples for tau-regression"))
        }
    	
        bar<- 1- K_beta/ K_gamma
    	
        return(bar)
    }
}
