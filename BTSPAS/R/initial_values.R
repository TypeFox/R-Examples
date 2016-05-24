## 2014-09-01 CJS bug in init.muLogTT which gives log(0) if a release has all recoveries only in initial strataum of
#                 of release. In those cases, I set the mean to those values that are not infinite
## 2013-12-18 CJS Any init.epsilon that correspond to logitP.fixed (typically to -10 or 0 on the p scale) must be set to NA
## 2013-09-04 CJS if any init.epsilon are NA, then set it to the mean of the non-missing values.
##                If any of n1, m2, u2 are missing set to average (but m2 <= n1). This tries to keep
##                OpenBugs from wandering off too far and generating nonsense values.
## 2012-02-01 CJS added na.rm=TRUE in computation of pScale to avoid passing NA
## 2011-05-15 CJS limited the etaU=log(U) to a maximum of 15 which corresponds to around 400,000,000 fish. 
## 2011-05-09 CJS subtle bug with initial values of epsilon where if fixed values for logitP at the end of the
##                experiment, then the initial values for epsilon must be truncated
## Function to generate initial values for each chain of the MCMC algorithm.
##
## Every model requires different initial values, though much of the
## code can be reused.

genInitsTTln <-
    function(n1,m2,u2){
        ## Generate initial parameters for log-normal travel time model
        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        m2dot1 <- apply(m2[,1:Nstrata.cap],1,sum)

        init.muLogTT <- rep(NA,Nstrata.rel)
        tmp1 <- (m2dot1>0)
#       init.muLogTT[tmp1] <- log((m2[tmp1,1:Nstrata.cap] %*% 1:Nstrata.cap)/(m2dot1[tmp1]) - (1:Nstrata.rel)[tmp1])
        init.muLogTT[tmp1] <- log(pmax(1, (m2[tmp1,1:Nstrata.cap] %*% 1:Nstrata.cap)/(m2dot1[tmp1]) - (1:Nstrata.rel)[tmp1]))  # 2014-09-01 added the pmax(1, xxx) to avoid taking log of zero

        init.muLogTT[!tmp1] <- mean(init.muLogTT[tmp1])

        init.xiMu <- mean(init.muLogTT)
        init.tauMu <- 1/max(0.2,var(init.muLogTT))  # avoid variances that are zero

        init.etasdLogTT <- log(rep(.5,Nstrata.rel))  # note that log (sd(log travel time)) is being modelled
        init.xiSd <- mean(init.etasdLogTT)
        init.tauSd <- 1

        return(list(muLogTT=init.muLogTT,
                    xiMu=init.xiMu,
                    tauMu=init.tauMu,
                    xiSd=init.xiSd,
                    tauSd=init.tauSd,
                    etasdLogTT=init.etasdLogTT))
    }

genInitsTTnp <-
    function(n1,m2,u2,Delta.max){
        ## Generate initial parameters for non-parametric travel time model

        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        ## Compute empirical theta matrix
        init.Theta <- t(sapply(1:Nstrata.rel,function(i){
            if(all(is.na(m2[i,])) || sum(m2[i,])==0)
                return(rep(NA,Delta.max+1))

            else{
                thetatmp <- pmax(.01,
                                 pmin(m2[i,-(Delta.max+2)]/sum(m2[i,-(Delta.max+2)],na.rm=TRUE),
                                      .99,na.rm=TRUE))  # CJS 2011-02-16
                return(thetatmp/sum(thetatmp))
            }
        }))

        ## Compute initial r
        init.delta <- t(apply(as.matrix(init.Theta[,-(Delta.max+1)]),1,   # CJS 2011-02-16 as.matrix added
                              function(theta){    # CJS fixed -(Delta.max+1)
                                  if(length(theta) == 1){theta}
                                  else {theta/(1-c(0,cumsum(theta[-Delta.max])))}
                              }))

        init.r <- log(init.delta)

        ## mean and standard deviation of transition probabilties
        init.muTT <- apply(logit(init.delta),2,mean,na.rm=TRUE)
        init.sdTT <- sd(as.vector(t(logit(init.delta)))-init.muTT,na.rm=TRUE)

        return(list(muTT=init.muTT,
                    tauTT=1/init.sdTT^2,
                    r=init.r,Theta=init.Theta))
    }

genInitValsChain <-
    function(model,
             n1,                          # Individuals marked per strata at first location
             m2,                          # Individuals recovered at second location
             u2,                          # (List of) unmarked individuals captured per strata with a single spline
             Delta.max=NULL,              # Max travel time for NP model
             logitP.cov,                  # Covariate matrix for capture probabilities
             logitP.fixed=NULL,           # Which logitP are fixed (typically to zero)?
             SplineDesign,                # (List of) design matrix(ces) for splines
             hatch.after=NULL,            # Data of release for hatchery fish in model with two splines
             pScale=1){

        #cat("\n** genInitValsChain \n")
        #browser()

        ## Generate initial values for a single chain
        Nstrata.rel <- length(n1)
        Nstrata.cap <- length(u2)

        inits <- list()                   # Create empty list of initial values

        ## 1) Travel time parameters (for non-diagonal models only)
        if(model %in% "TSPNDE"){
            inits <- append(inits,genInitsTTln(n1,m2,u2))
        }
        if(model %in% "TSPNDENP"){
            inits <- append(inits,genInitsTTnp(n1,m2,u2,Delta.max))
        }

        ## 2) Capture probabilities
        ## 2.1) Compute initial logit capture probabilities
        if(model %in% c("TSPDE","TSPDE-WHchinook","TSPDE-WHsteel")){
            init.P <- (m2+1)/(n1+2) * pScale
        }
        else if(model %in% c("TSPNDE")){
            ## Compute expected number of marked fish in each cell
            Theta <- t(sapply(1:Nstrata.rel,function(i){
                tmp <- pnorm(log(i:Nstrata.cap),inits$muLogTT[i],exp(inits$etasdLogTT[i]))
                c(rep(0,(i-1)),tmp - c(0,tmp[-(Nstrata.cap - (i-1))]))
            }))

            M <- Theta * n1

            m2dot2 <- apply(m2[,1:Nstrata.cap],2,sum)
            init.P <- (m2dot2 + 1)/(apply(M,2,sum) + m2dot2 + 1) * pScale
        }
        else if(model %in% c("TSPNDENP")){
            ## Compute expected number of marked fish in each cell
             N2 <- lapply(1:Nstrata.rel,function(i) inits$Theta[i,]*n1[i])

            ## Compute expected number of marked fish in each capture strata
            n2 <- sapply(1:Nstrata.cap,function(i){
                n2tmp <- 0
                for(j in max(i-Delta.max,1):min(i,Nstrata.rel))
                    n2tmp <- N2[[j]][i-j+1] + n2tmp
                n2tmp
            })

            m2dot2 <- sapply(1:Nstrata.cap,function(i){
                m2tmp <- 0
                for(j in max(i-Delta.max,1):min(i,Nstrata.rel))
                    m2tmp <- m2[j,i-j+1] + m2tmp
                m2tmp
            })

            init.P <- (m2dot2 + 1)/(n2 + m2dot2 + 1) * pScale

        }

        init.P <- pmax(.00001,pmin(init.P,.99999))  # constrain p to the interval (.00001, .99999)
        init.P[is.na(init.P)] <- mean(init.P,na.rm=TRUE) # remove missing values
        init.P[!is.na(logitP.fixed)] <- NA # remove fixed values from initial vector
        init.logitP <- logit(init.P)    # Compute the logit

        ## 2.2) Compute associated coefficients for design matrix
        init.beta.logitP <- as.vector(lm(init.logitP ~ logitP.cov - 1)$coeff)

        ## 2.3) Set variance for hierarchical model of capture probabilities
        if(length(init.beta.logitP)==1)
            init.tauP <- 1/var(init.logitP - logitP.cov*init.beta.logitP,na.rm=TRUE)
        else
            init.tauP <- 1/var(init.logitP - logitP.cov %*% init.beta.logitP,na.rm=TRUE)

        init.beta.logitP <- c(init.beta.logitP, 0)   # add one extra element so that single beta is still written as a vector

        init.beta.logitP[is.na(init.beta.logitP)] <- 0

        inits <- append(inits,list(beta.logitP=init.beta.logitP,tauP=as.numeric(init.tauP)))

        ## 3) Numbers of unmarked individuals per strata (where u2 observed)
        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            init.U <- ceiling((u2+1)/init.P)
        }

        ## Option 2: Chinook model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHchinook")){
            init.U.W <- ceiling((u2$W+1)/init.P)
            init.U.H <- ceiling((u2$H+1)/init.P)
            init.U.H[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery
        }

        ## Option 3: Steelhead model with separate splines for wild, wild YOY, and hatchery fish
        if(model %in% c("TSPDE-WHsteel")){
            init.U.W.YoY <- ceiling((u2$W.YoY+1)/init.P)
            init.U.W.1 <- ceiling((u2$W.1+1)/init.P)
            init.U.H.1 <- ceiling((u2$H.1+1)/init.P)
            init.U.H.1[1:hatch.after] <- 0   # no hatchery fish prior to release from hatchery
        }

        ## 4) Spline coefficients

        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            ## 4.1) Fit Spline to strata with u2 observed
            tmp1 <- !is.na(init.U)
            init.bU <- lm(log(init.U[tmp1]) ~ SplineDesign[tmp1,]-1)$coeff
            init.bU[is.na(init.bU)] <- mean(init.bU,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.2) Compute variance of second differences between coefficients
            tmp2 <- 3:length(init.bU)
            sigmaU <- sd(init.bU[tmp2]-2*init.bU[tmp2-1]+init.bU[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU=init.bU,tauU=init.tauU))
        }

        ## Option 2: Chinook model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHchinook")){
            ## 4.1.a) Fit spline to wild fish
            tmp1.W <- !is.na(init.U.W)
            init.bU.W <- lm(log(init.U.W[tmp1.W]) ~ SplineDesign$W[tmp1.W,]-1)$coeff
            init.bU.W[is.na(init.bU.W)] <- mean(init.bU.W,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.b) Fit spline to hatchery fish
            tmp1.H <- c(rep(FALSE,hatch.after),!is.na(init.U.H[-(1:hatch.after)]))
            init.bU.H <- lm(log(init.U.H[tmp1.H]) ~ SplineDesign$H[tmp1.H,]-1)$coeff
            init.bU.H[is.na(init.bU.H)] <- mean(init.bU.H,na.rm=TRUE)       # Fix any coefficients that can't be c

            ## 4.2) Variance of second differences between coefficients (use only wild fish to initialize)
            tmp2 <- 3:length(init.bU.W)
            sigmaU <- sd(init.bU.W[tmp2]-2*init.bU.W[tmp2-1]+init.bU.W[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU.W=init.bU.W,bU.H=init.bU.H,tauU=init.tauU))
        }

        ## Option 3: Steelhead model with separate splines for wild and hatchery fish
        if(model %in% c("TSPDE-WHsteel")){
            ## 4.1.a) Fit spline to wild YoY fish
            tmp1.W.YoY <- !is.na(init.U.W.YoY)
            init.bU.W.YoY <- lm(log(init.U.W.YoY[tmp1.W.YoY]) ~ SplineDesign$W.YoY[tmp1.W.YoY,]-1)$coeff
            init.bU.W.YoY[is.na(init.bU.W.YoY)] <- mean(init.bU.W.YoY,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.b) Fit spline to wild 1 fish
            tmp1.W.1 <- !is.na(init.U.W.1)
            init.bU.W.1 <- lm(log(init.U.W.1[tmp1.W.1]) ~ SplineDesign$W.1[tmp1.W.1,]-1)$coeff
            init.bU.W.1[is.na(init.bU.W.1)] <- mean(init.bU.W.1,na.rm=TRUE)       # Fix any coefficients that can't be computed

            ## 4.1.c) Fit spline to hatchery fish
            tmp1.H.1 <- c(rep(FALSE,hatch.after),!is.na(init.U.H.1[-(1:hatch.after)]))
            init.bU.H.1 <- lm(log(init.U.H.1[tmp1.H.1]) ~ SplineDesign$H.1[tmp1.H.1,]-1)$coeff
            init.bU.H.1[is.na(init.bU.H.1)] <- mean(init.bU.H.1,na.rm=TRUE)       # Fix any coefficients that can't be c

            ## 4.2) Variance of second differences between coefficients (use only wild YoY fish to initialize)
            tmp2 <- 3:length(init.bU.W.YoY)
            sigmaU <- sd(init.bU.W.YoY[tmp2]-2*init.bU.W.YoY[tmp2-1]+init.bU.W.YoY[tmp2-2])
            init.tauU <- 1/sigmaU^2

            inits <- append(inits,list(bU.W.YoY=init.bU.W.YoY,bU.W.1=init.bU.W.1,
                                       bU.H.1=init.bU.H.1,tauU=init.tauU))
        }

        ## 5) Variance about spline
        ## Option 1: Models with one spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            sigmaeU <- sd(log(init.U+1) - SplineDesign %*% init.bU,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        ## Option 2: Chinook models with two splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHchinook")){
            sigmaeU <- sd(log(init.U.W+1) - SplineDesign$W %*% init.bU.W,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        ## Option 3: Steelhead models with three splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHsteel")){
            sigmaeU <- sd(log(init.U.W.YoY+1) - SplineDesign$W.YoY %*% init.bU.W.YoY,na.rm=TRUE)
            init.taueU <- 1/sigmaeU^2
        }

        inits <- append(inits,list(taueU=init.taueU))

        ## 6) Initialize missing U values by fitting spline and generating errors
        ## Option 1: Models with only 1 spline
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            if(sum(!tmp1)>0)
                init.U[!tmp1] <- ceiling(exp(as.vector(SplineDesign[!tmp1,] %*% init.bU)
                                             + rnorm(sum(!tmp1),0,sigmaeU))) + 1
            init.etaU <- pmin(log(init.U),20)  # limit the initial values to reasonable values

            inits <- append(inits,list(etaU=init.etaU))
        }

        ## Option 2: Chinook models with two splines
        if(model %in% c("TSPDE-WHchinook")){
            ## Wild fish
            if(sum(!tmp1.W)>0)
                init.U[!tmp1.W] <- ceiling(exp(as.vector(SplineDesign$W[!tmp1.W,] %*% init.bU.W)
                                               + rnorm(sum(!tmp1.W),0,sigmaeU))) + 1
            init.etaU.W <- pmin(log(init.U.W), 15) # limit the initial values to reasonable values

            ## Hatchery fish
            tmp2.H <- tmp1.H[-(1:hatch.after)]
            if(sum(!tmp2.H)>0)
                init.U[!tmp2.H] <- ceiling(exp(as.vector(SplineDesign$H[!tmp2.H,] %*% init.bU.H)
                                               + rnorm(sum(!tmp2.H),0,sigmaeU))) + 1

            init.etaU.H <- c(rep(NA,hatch.after),pmin(20,log(init.U.H[tmp1.H])))
      

            inits <- append(inits,list(etaU.W=init.etaU.W,etaU.H=init.etaU.H))
        }

        ## Option 3: Steelhead models with three splines
        if(model %in% c("TSPDE-WHsteel")){
            ## Wild YoY fish
            if(sum(!tmp1.W.YoY)>0)
                init.U[!tmp1.W.YoY] <- ceiling(exp(as.vector(SplineDesign$W.YoY[!tmp1.W.YoY,] %*% init.bU.W.YoY)
                                                   + rnorm(sum(!tmp1.W.YoY),0,sigmaeU))) + 1
            init.etaU.W.YoY <- pmin(20,log(init.U.W.YoY))


            ## Wild 1 fish
            if(sum(!tmp1.W.1)>0)
                init.U[!tmp1.W.1] <- ceiling(exp(as.vector(SplineDesign$W.1[!tmp1.W.1,] %*% init.bU.W.1)
                                                 + rnorm(sum(!tmp1.W.1),0,sigmaeU))) + 1
            init.etaU.W.1 <- pmin(20,log(init.U.W.1))

            ## Hatchery 1 fish
            tmp2.H.1 <- tmp1.H.1[-(1:hatch.after)]
            if(sum(!tmp2.H.1)>0)
                init.U[!tmp2.H.1] <- ceiling(exp(as.vector(SplineDesign$H.1[!tmp2.H.1,] %*% init.bU.H.1)
                                                 + rnorm(sum(!tmp2.H.1),0,sigmaeU))) + 1

            init.etaU.H.1 <- c(rep(NA,hatch.after),pmin(20,log(init.U.H.1[tmp1.H.1])))

            inits <- append(inits,list(etaU.W.YoY=init.etaU.W.YoY,
                                       etaU.W.1=init.etaU.W.1,etaU.H.1=init.etaU.H.1))
        }

        ## 7) Transform initial values for logitP to initial values for epsilon
        #     If some of the logitP are fixed, you need to set the corresponding value of epsilon to NA
        #     This is done at the end of these possible model choices.
        ## Option 1: Models with only one spline
        #cat("GenInitVals \n")
        #browser()
        if(model %in% c("TSPDE","TSPNDE","TSPNDENP")){
            #cat("GenInitVals - setting epsilon: ", model, "\n")
            #browser()
            init.epsilon <- init.logitP - log(u2 + 1) + inits$etaU
            # subtle problem. If the logitP.fixed include elements at the end of the experiment
            # then init.epsion needs to be truncated at the end, otherwise OPENBugs gets upset
            if(length(logitP.fixed)>0){
               for(i in length(logitP.fixed):1){
                  if( is.na(logitP.fixed[i])){break}
                  init.epsilon <- init.epsilon[-length(init.epsilon)] # drop last term
               }
            }
        }

        ## Option 2: Chinook models with two splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHchinook")){
            init.epsilon <- init.logitP - log(u2$W + 1) + inits$etaU.W
        }

        ## Option 3: Steelhead models with three splines -- use only wild fish to initialize
        if(model %in% c("TSPDE-WHsteel")){
            init.epsilon <- init.logitP - log(u2$W.YoY + 1) + inits$etaU.W.YoY
        }

        ## Change any missing epsilon values to the mean of the epsilon unless these were from 
        ## logitP values that were fixed. In those cases, the epsilon must remain as missing
        init.epsilon[is.na(init.epsilon)] <- mean(init.epsilon, na.rm=TRUE)  ## CJS 2013-09-04
        init.epsilon[!is.na(logitP.fixed)] <- NA  ## CJS 2013-12-17 (if logitP is fixed, don't initialize epsilon
        inits <- append(inits,list(epsilon=c(init.epsilon)))

        ## Remove working objects from the initial values
        if(model=="TSPNDENP"){
            inits$Theta <- NULL
        }

        ##9. Generate initial values for missing n1, m2, or u2 as the average of the other values (rounded to integers)
        if(model %in% c("TSPDE","TSPDE-WHchinook","TSPNDE","TSPNDENP","TSPDE-WHsteel")){
           init.n1 <- rep(NA, length(n1))
           init.n1[is.na(n1)] <- round(mean(n1, na.rm=TRUE))
           if(any(is.na(n1))){inits <- append(inits, list(n1=init.n1))}
        }
 
        if(model %in% c("TSPDE","TSPDE-WHchinook","TSPDE-WHsteel")){ 
           init.m2 <- rep(NA, length(m2))
           init.m2[is.na(m2)] <- round(pmin(n1[is.na(m2)],mean(m2, na.rm=TRUE)))
           if(any(is.na(m2))){inits <- append(inits, list(m2=init.m2))}
        }
        if(model %in% c("TSPNDE","TSPNDENP")){ 
           # not sure how to initialize bad m2 values for the non-diagonal case
           # because m2 is a full matrix with elements arranges diagonally
        }

        if(model %in% c("TSPDE", "TSPNDE","TSPNDENP")){
           init.u2 <- rep(NA, length(u2))
           init.u2[is.na(u2)] <- round(mean(u2, na.rm=TRUE))
           if(any(is.na(u2))){inits <- append(inits, list(u2=init.u2))}
        }
        if(model %in% c("TSPDE-WHchinook")){  # u2 is a list with components W and H, A and N
           init.u2.A <- rep(NA, length(u2$A))
           init.u2.A[is.na(u2$A)] <- pmin(init.U.H[is.na(u2$A)], round(median(u2$A, na.rm=TRUE)))
           if(any(is.na(u2$A))){inits <- append(inits, list(u2.A=init.u2.A))}
           init.u2.N <- rep(NA, length(u2$N))
           init.u2.N[is.na(u2$N)] <- pmin(init.U.W[is.na(u2$N)],round(median(u2$N, na.rm=TRUE))) # This is too strict as some hatchery have no clips
           if(any(is.na(u2$N))){inits <- append(inits, list(u2.N=init.u2.N))}
        }
        if(model %in% c("TSPDE-WHsteel")){  # u2 is a list with components W.YoY, W.1, H.1
           init.u2.W.1 <- rep(NA, length(u2$W.1))
           init.u2.W.1[is.na(u2$W.1)] <- pmin(init.U.W.1[is.na(u2$W.1)], round(median(u2$W.1, na.rm=TRUE)))
           if(any(is.na(u2$W.1))){inits <- append(inits, list(u2.W.1=init.u2.W.1))}
           init.u2.W.YoY <- rep(NA, length(u2$W.YoY))
           init.u2.W.YoY[is.na(u2$W.YoY)] <- pmin(init.U.W.YoY[is.na(u2$W.YoY)], round(median(u2$W.YoY, na.rm=TRUE)))
           if(any(is.na(u2$W.YoY))){inits <- append(inits, list(u2.W.YoY=init.u2.W.YoY))}
           init.u2.H.1 <- rep(NA, length(u2$H.1))
           init.u2.H.1[is.na(u2$H.1)] <- pmin(init.U.H.1[is.na(u2$H.1)], round(median(u2$H.1, na.rm=TRUE)))
           if(any(is.na(u2$H.1))){inits <- append(inits, list(u2.H.1=init.u2.H.1))}
        }

        return(inits)

    }

genInitVals <-
    function(model,
             n1,                          # Individuals marked per strata at first location
             m2,                          # Individuals recovered at second location
             u2=NULL,                     # (List of) unmarked individuals captured
             Delta.max=NULL,              # Max travel time for NP model
             logitP.cov,                  # Covariate matrix for capture probabilities
             logitP.fixed=NULL,           # Which values of logitP are fixed (typically at zero)?
             SplineDesign,                # (List of) desgin matrices for spline for models
             hatch.after=NULL,            # Data of release for hatchery fish in model with two splines
             n.chains=3,
             pStep=5){                   # Relative change in p between chains


      ## Determine maximum scaling factor in order to avoid tauP=Inf
      pScaleMax <- 1/(1.1*sum(m2,na.rm=TRUE)/sum(n1,na.rm=TRUE))
      
      ## Generate initial values for n.chains
      inits <- lapply(1:n.chains,function(i){
        ## Compute scaling factor for ith chain
        #cat("\n*** gen init values ***\n")
        #browser()
        pScale <- min(pStep ^(-(n.chains-1)/2 + (i-1)),pScaleMax,1,na.rm=TRUE)

        ## Generate initial values
        genInitValsChain(model,
                         n1,m2,u2,Delta.max,
                         logitP.cov,logitP.fixed,
                         SplineDesign,hatch.after,
                         pScale)
      })
      
      return(inits)
    }
