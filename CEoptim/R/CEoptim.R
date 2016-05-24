CEoptim <- function(f, f.arg=NULL, maximize=FALSE, continuous=NULL, discrete=NULL, N=100L, rho=0.1, iterThr=1e4L, noImproveThr=5, verbose=FALSE ) {

### Parse and check arguments.

    ## f should be a function
    if (!is.function(f)) stop("Argument 1 (f) should be a function.")

    ## Argument <continuous> should be a list specifying the sampling
    ## distribution of the continuous optimization
    ## variables. Incidentally it will define the number of such
    ## variables.

    ## The continuous sampling parameters and their defaults.
    
    ctsDefaults<- list(mean=NULL, sd=NULL, conMat=NULL, conVec=NULL,smoothMean=1, smoothSd=1, sdThr=0.001)

    if (!is.null(continuous)) {
        if (!is.list(continuous))
            stop('Argument 2 (continuous) should be a list specifying the sampling distribution of the continuous optimization variables.')
        nameCon <- names(ctsDefaults)   #name inner list
        nameInput <- names(continuous)
        if (length(noNms <- nameInput[!nameInput %in% nameCon])!=0) {
            stop("unknown parameter(s) in continuous: ",
                 paste(noNms, collapse = ", "))
        }  else ctsDefaults[nameInput]<-continuous
    }

    mu0 <- ctsDefaults$mean
    sigma0 <- ctsDefaults$sd
    A <- ctsDefaults$conMat
    b <- ctsDefaults$conVec
    alpha <- ctsDefaults$smoothMean
    beta <- ctsDefaults$smoothSd
    eps <- ctsDefaults$sdThr
    
    ## Argument <discrete> should be a list specifying the sampling
    ## distribution of the discrete optimization variables. Incidentally
    ## it will define the number of such variables.

    ## The discrete sampling parameters and their defaults.
    
    discDefaults<-list(categories=NULL,probs=NULL, smoothProb=1,probThr=0.001)

    if (!is.null(discrete)) {
        if (!is.list(discrete))
            stop('Argument 3 (discrete) should be a list specifying the sampling distribution of the discrete optimization variables.')
        nameDis <-names(discDefaults)
        nameDinput <-names(discrete)
        if (length(noNms <- nameDinput[!nameDinput %in% nameDis])!=0) {
            stop("unknown parameter(s) in discrete: ", paste(noNms, collapse = ", "))
        } else discDefaults[nameDinput]<-discrete
    }

    categories<-discDefaults$categories
    tau0<- discDefaults$probs
    gamma <- discDefaults$smoothProb
    eta <- discDefaults$probThr


    r <- NULL
    unspecMu <- NULL
    if (!is.null(mu0) || !is.null(sigma0)) {
        if (is.null(mu0) || is.null(sigma0))
            stop("If mu0 is specified so must be sigma0, and vice versa.")
                                        #    if (!is.matrix(sigma0))
        Sigma0 <- diag(sigma0,nrow=length(sigma0))^2 #QB Sigma0
        if (length(mu0)!=length(sigma0))
            stop("Arguments \"mu\" and  \"sigma\" must be of same length.")
        else r <- length(mu0)
        ## The user might give NA as distribution parameter for some variables.
        unspecMu <- is.na(mu0)
    }

    ## If A is specified it should be a p-column matrix.
    if (!is.null(A)) {
        if (!is.matrix(A)) stop("Argument \"A\" should be a matrix.")
        p <- dim(A)[2]          # problem cts dimension
        c <- dim(A)[1]          # no. cts constraints
                                        # Test conformity with mu0 and sigma0.
        if (!is.null(r) && r!=p)
            stop('Number of continuous variables implied by "A" argument does not match that implied by mu0 and sigma0.')
                                        # Set unspecMu if necessary.
        if (is.null(unspecMu)) unspecMu <- rep(TRUE,p)
                                        # Test b given
        if (is.null(b)) stop('Argument "A" requires argument "b"')
        if (length(b)!=c) stop('Argument "b" not same length as no. rows A')
    }
    else {
        p <- r
        if (is.null(p)) p <- 0L
        if (!is.null(b)) stop('Argument "b" requires argument "A"')
    }

    ## If there are continuous vars then they must have an initial
    ## distribution of some sort.
    ## Could test whether the constraints define a bounded search region.
    if (any(unspecMu))
        stop("Mu and Sigma must be specified for continuous variables.")

    ## if categories is specified it should be a vector of integers

    if (!is.null(categories) && is.null(tau0)) {
        if (!is.vector(categories) || !identical(typeof(categories),'integer'))
            stop("Argument \"categories\" should be an integer vector.")
        q <- length(categories)
        tau0 <- list();
        for (c in categories) {
            tau0 <- c(tau0,list(rep(1.0/c,c)))
        }
    }

    if (!is.null(tau0)){
        if(!is.list(tau0))
            stop("\"tau\" in Argument discrete should be a list")
        if(is.null(categories))
            categories<-mapply(length,tau0)
    }
    q <- length(categories)
    ## There should be at least one argument to f.
    if (p+q==0L)
        stop('There should be at least one argument discrete or continuous to f.')
    ## Elite prop. rho should be between 0 and 1.
    if (!is.numeric(rho) || !is.vector(rho) || length(rho)!=1
        || rho>=1.0 || rho<=0.0)
        stop("Argument \"rho\" should be a number between 0 and 1\n")
    nElite <- round(N*rho)

    ## Update smoothing pars. alpha, beta and gamma should be between 0 and 1.
    if (any(!is.numeric(c(alpha,beta,gamma)))
        || any(c(alpha,beta,gamma)<0) || any(c(alpha,beta,gamma)>1))
        stop('Arguments alpha, beta and gamma should be between 0 and 1.')

    ## Echo args for debugging.

    if (verbose) {
        cat('Number of continuous variables:',p,' \n')
        cat('Number of discrete variables:',q,'\n')
        cat('conMat=','\n')
        print(A)
        cat('conVec=','\n')
        print(b)
        cat('smoothMean:',alpha,'smoothSd:',beta,'smoothProb:',gamma,'\n')
        cat('N:',N,'rho:',rho,'iterThr:',iterThr,'sdThr:',eps,'probThr',eta,'\n')
    }
    ## Create a wrapper for f that takes two arguments, one for the
    ## continuous optimization variables, one for the discrete.
    s = (-1)^maximize
    if (p==0L) ff <- function(XC,XD) { return(s*do.call(f,
            c(list(XD),f.arg))) }
    else if (q==0L)  ff <- function(XC,XD) { s*do.call(f,c(list(XC),f.arg)) }
    else  ff <- function(XC,XD) { s*do.call(f,c(list(XC,XD),f.arg)) }

    ## Generate an initial sample.
    Xc <- matrix(nrow=N,ncol=p)         #cts portion of sample
    if (p>0) {
        Xc <- rtmvnorm(N,mu0,Sigma0,A,b)$X  #QB Sigma0 Variance
    }

    Xd <- matrix(nrow=N,ncol=q)

    if (q>0) {

        for (i in 1:q) {
            Xd[,i] <- sample(0:(categories[i]-1),N,replace=T,prob=tau0[[i]])
        }
    }


    ## Evaluate objective function over initial sample
    Y <- mapply(ff,lapply(1:N,function(j)Xc[j, ,drop=F]),
                lapply(1:N,function(k)Xd[k, ,drop=F]))

    ## Identify elite.
    IX <- sort(Y,index.return=TRUE,decreasing=F)$ix
    elite <- IX[1:nElite]
    

    ## Estimate sampling distributions.

    if (p>0) {
        ## Estimate a normal distribution with smoothing
        mu <- colMeans(Xc[elite, ,drop=F])*alpha + mu0*(1-alpha)
        sigma <- apply(Xc[elite, ,drop=F],MARGIN=2,FUN=sd)*beta + sigma0*(1-beta) #QB Smoothing parameter for the standard deviations 
        Sigma <- diag(sigma,nrow=p)^2
    }

    if (q>0) {
        
        counts <- lapply(split(Xd[elite, ,drop=F],col(Xd[elite, ,drop=F])),
                         table)
       
        tau <- list()
        for (i in 1:q) {
            v <- rep(0,categories[i])
            v[1+as.numeric(dimnames(counts[[i]])[[1]])] <-
                as.vector(counts[[i]])
            tau[[i]] <- v/nElite
        }

        ## Combine tau and tau0.
        tau <- mapply('+',lapply(tau,'*',gamma),lapply(tau0,'*',1-gamma),
                      SIMPLIFY=F)
    }


    iter <- 0
    ctsOpt <- Xc[elite[1],]
    disOpt <- Xd[elite[1],]
    optimum <- Y[elite[1]]
    gammat <-Y[elite[nElite]]
    ceprocess <- NULL
    diffopt <- Inf
    CEstates<- NULL
    probst<-list()
### Main loop -- test termination conditions
    while (iter < iterThr && diffopt!=0 &&
           (( p>0 && max(sigma) > eps) ||
            (q>0 && max(1.0-sapply(tau,max)) > eta))) {
     
      
      
      CEt<- NULL
      CEt<- c(iter,optimum*s,gammat*s)
      if(p>0){
        CEt<- c(CEt,mu,max(sigma))}
      if(q>0){
        CEt<- c(CEt,max(1.0-sapply(tau,max)))       
        namet<-paste("probs",iter,sep="")
        assign(namet,tau)
        probst[[iter+1]]<-get(namet)}
      
      CEstates<- rbind(CEstates,CEt)
      
      if (verbose) {  
        cat('iter:',iter,
            ' opt:',optimum*s)
       
        if(p>0){
          cat(' maxSd:',max(sigma))}
        
        if(q>0){
          cat(' maxProbs:', max(1.0-sapply(tau,max)))          
        }
        cat('\n')
        
      }
      
      
#       if (verbose) {
#             CEt<- NULL
#             
#             cat('iter:',iter,
#                 ' opt:',optimum*s)
#             CEt<- c(iter,optimum*s,gammat*s)
#             if(p>0){
#                 cat(' maxSd:',max(sigma))
#                 CEt<- c(CEt,mu,max(sigma))}
#             
#             if(q>0){
#                 cat(' maxProbs:', max(1.0-sapply(tau,max)))
#                 CEt<- c(CEt,max(1.0-sapply(tau,max)))
#                 
#                 namet<-paste("probs",iter,sep="")
#                 assign(namet,tau)
#                probst[[iter+1]]<-get(namet)
#                
#             }
#                      
#             CEstates<- rbind(CEstates,CEt)
#             cat('\n')
# 
#         }
#          
       
        ## Generate sample and evaluate objective function

        if (p>0) {
            ## Sample truncated normal distributions.
            Xc <- rtmvnorm(N,mu,Sigma,A,b)$X  #QB change sigma to Sigma
        }
        if (q>0) {
            ## sample categorical distributions.
            ## tjb should use some sort of "apply" instead of for loop.
            for (i in 1:q) {
                Xd[,i] <- sample(0:(categories[i]-1),N,replace=T,prob=tau[[i]])
            }
        }

        Y <- mapply(ff,lapply(1:N,function(j)Xc[j, ,drop=F]),
                    lapply(1:N,function(k)Xd[k, ,drop=F]))

        ## Identify elite.
        IX <- sort(Y,index.return=TRUE,decreasing=F)$ix
        elite <- IX[1:nElite]

        ## test for new optimum
        if (Y[elite[1]] < optimum ) {
           
            if (p>0) ctsOpt <- Xc[elite[1],]
            if (q>0) disOpt <- Xd[elite[1],]
            optimum <- Y[elite[1]]
            if(Y[elite[nElite]]<gammat)
              gammat<-Y[elite[nElite]]
        }# change here

        
        #ce process
        ceprocess <- c(ceprocess,optimum)
        if (iter > noImproveThr ){
            diffopt<-sum(abs(ceprocess[(iter-noImproveThr):(iter-1)]-optimum))
        }


        ## Reestimate sampling distributions.
        if (p>0) {
            mu <- colMeans(Xc[elite, ,drop=F])*alpha + mu*(1-alpha)
            
            sigma <- apply(Xc[elite, ,drop=F],MARGIN=2,FUN=sd)*beta + sigma*(1-beta)
            Sigma <- diag(sigma,nrow=p)^2

        }

        if (q>0) {
            ## Reestimate multivariate categorical distribution with smoothing.
            tau0 <- tau
            counts <- lapply(split(Xd[elite, ,drop=F],col(Xd[elite, ,drop=F])),
                             table)
            tau <- list()
            for (i in 1:q) {
                v <- rep(0,categories[i])
                v[1+as.numeric(dimnames(counts[[i]])[[1]])] <-
                    as.vector(counts[[i]])
                tau[[i]] <- v/nElite
            }
            ## Combine tau and tau0.
            tau <- mapply('+',lapply(tau, '*',gamma),lapply(tau0,'*',1-gamma),
                          SIMPLIFY=F)
        }
        iter <- iter+1
    }
      
     

     if(iter==iterThr)
       convergence="Not converged" 
     else if(diffopt==0)
       convergence=paste("Optimum did not change for",noImproveThr,"iterations")
     else
       convergence="Variance converged"
    if(verbose){
    rownames(CEstates)<-c()
    if(p>0&&q==0)
        colnames(CEstates)<- c("iter","optimum","gammat",paste("mean",1:p,sep=""),"maxSd")
    else if(p==0 && q>0)
        colnames(CEstates)<- c("iter","optimum","gammat","maxProbs")
    else if(p>0 && q>0)
      colnames(CEstates)<- c("iter","optimum","gammat",paste("mean",1:p,sep=""),"maxSd","maxProbs")
    }

   out<-list(optimizer=list(continuous=ctsOpt,discrete=disOpt),
          optimum=s*optimum,termination=list(niter=iter, nfe=iter*N, convergence=convergence),
          states=CEstates,states.probs=probst)
   class(out)<- "CEoptim"
   out

}

