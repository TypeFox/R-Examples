initiate.startValues <- function(Y, lin.pred, data, model, cluster=NULL,
                                 beta1=NULL, beta2=NULL, beta3=NULL, beta=NULL,
                                 gamma.ji=NULL,
                                 theta=NULL,
                                 V.j1=NULL, V.j2=NULL, V.j3=NULL, V.j = NULL,
                                 WB.alpha=NULL, WB.kappa=NULL,
                                 PEM.lambda1=NULL, PEM.lambda2=NULL, PEM.lambda3=NULL, PEM.lambda=NULL,
                                 PEM.s1=NULL, PEM.s2=NULL, PEM.s3=NULL, PEM.s=NULL,
                                 PEM.mu_lam=NULL, PEM.sigSq_lam=NULL,
                                 MVN.SigmaV=NULL, Normal.zeta=NULL,
                                 DPM.class=NULL, DPM.tau=NULL)
{
    
    
    
    
    ## BayesSurv
    
    if(length(model)==1)
    {
        print(paste("Start values are initiated for univariate ", model, " model...", sep = ""), cat("\n"))
        ##
        if(!is.null(cluster)) print(paste("Warning: 'cluster' is not required for ", model, " model so it is ignored", sep = ""))
        
        ##
        n <- nrow(Y)
        p <- ncol(model.frame(lin.pred, data=data))
        
        
        ##
        if(!is.null(beta1)) stop(paste("'beta1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(beta2)) stop(paste("'beta2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(beta3)) stop(paste("'beta3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(gamma.ji)) stop(paste("'gamma.ji' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(theta)) stop(paste("'theta' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class)) stop(paste("'DPM.class' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.tau)) stop(paste("'DPM.tau' is for univariate models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda1)) stop(paste("'PEM.lambda1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda2)) stop(paste("'PEM.lambda2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda3)) stop(paste("'PEM.lambda3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s1)) stop(paste("'PEM.s1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s2)) stop(paste("'PEM.s2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s3)) stop(paste("'PEM.s3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        
        
        ##
        if(!is.null(beta)){
            if(length(beta) != p) stop(paste("Length of starting value for beta must be", p))
        }
        if(is.null(beta)) beta <- runif(p, -0.1, 0.1)
        
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 1) stop("Length of starting value for WB.alpha must be 1")
        }
        if(is.null(WB.alpha)) WB.alpha <- 1
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 1) stop("Length of starting value for WB.kappa should be 1")
        }
        if(is.null(WB.kappa)) WB.kappa <- 0.01
        
        ##
        if(!is.null(PEM.s))
        {
        	if(PEM.s[1]== 0) stop("'min(PEM.s)' must be greater than 0")
        } 
        if((!is.null(PEM.lambda)) + (!is.null(PEM.s)) == 1) stop("Both 'PEM.s' and 'PEM.lambda' must be specified as numeric vectors or as NULL")
        if(!is.null(PEM.lambda) & !is.null(PEM.s))
        {
            if(length(PEM.lambda) != length(PEM.s)) stop("Length of starting value for 'PEM.lambda' and 'PEM.s' should be equal")
        }
        if(is.null(PEM.lambda) & is.null(PEM.s))
        {
            uniqTime <- sort(unique(Y[Y[,2]==1,1]))
            Tmax <- max(uniqTime)
            if(length(uniqTime)>20)
            {
                PEM.s <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
            }else
            {
                PEM.s <- uniqTime
            }
            PEM.lambda <- runif(length(PEM.s), -3, -2)
        }
        ##
        if(!is.null(PEM.mu_lam)){
            if(length(PEM.mu_lam) != 1) stop("Length of starting value for PEM.mu_lam should be 1")
        }
        if(is.null(PEM.mu_lam)) PEM.mu_lam <- mean(PEM.lambda)
        ##
        if(!is.null(PEM.sigSq_lam)){
            if(length(PEM.sigSq_lam) != 1) stop("Length of starting value for PEM.sigSq_lam should be 1")
        }
        if(is.null(PEM.sigSq_lam)) PEM.sigSq_lam <- ifelse(length(PEM.lambda)==1, 0.1, var(PEM.lambda))
        
        ##
        start.common <- list(beta=beta)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        start.PEM   <- list(PEM.lambda=PEM.lambda, PEM.s=PEM.s, PEM.mu_lam=PEM.mu_lam, PEM.sigSq_lam=PEM.sigSq_lam, K=length(PEM.s)-1)
        
        ##
        if(model == "Weibull"){
            value <- list(common=start.common, WB=start.WB)
        }
        if(model == "PEM"){
            value <- list(common=start.common, PEM=start.PEM)
        }
        
    }
    
    
    
    
    ## BayesSurvcor
    
    if(length(model)==2 & class(lin.pred)=="formula")
    {
      print(paste("Start values are initiated for univariate ", model[1],"-", model[2], " model...", sep = ""), cat("\n"))
      ##
      if(is.null(cluster)) stop(paste("'cluster' must be given for ", model[1],"-", model[2], " model.", sep = ""))
      
        ##
        n <- nrow(Y)
        J <- length(unique(cluster))
        p <- ncol(model.frame(lin.pred, data=data))
        
        ##
        if(!is.null(beta1)) stop(paste("'beta1' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(beta2)) stop(paste("'beta2' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(beta3)) stop(paste("'beta3' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(gamma.ji)) stop(paste("'gamma.ji' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(theta)) stop(paste("'theta' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for semi-competing risks models so it must be specified as NULL for ", model[1],"-", model[2], " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda1)) stop(paste("'PEM.lambda1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda2)) stop(paste("'PEM.lambda2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.lambda3)) stop(paste("'PEM.lambda3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s1)) stop(paste("'PEM.s1' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s2)) stop(paste("'PEM.s2' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(PEM.s3)) stop(paste("'PEM.s3' is for semi-competing risks models so it must be specified as NULL for ", model, " model.", sep = ""))
        
        
        ##
        if(!is.null(beta)){
            if(length(beta) != p) stop(paste("Length of starting value for beta must be", p))
        }
        if(is.null(beta)) beta <- runif(p, -0.1, 0.1)
        
        ##
        if(!is.null(V.j)){
            if(length(V.j) != J) stop(paste("Length of starting values for V.j must be", J))
        }
        if(is.null(V.j)) V.j <- runif(J, -0.1, 0.1)
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 1) stop("Length of starting value for WB.alpha must be 1")
        }
        if(is.null(WB.alpha)) WB.alpha <- 1
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 1) stop("Length of starting value for WB.kappa should be 1")
        }
        if(is.null(WB.kappa)) WB.kappa <- 0.01
        
        ##
        if(!is.null(PEM.s))
        {
        	if(PEM.s[1]== 0) stop("'min(PEM.s)' must be greater than 0")
        }
        if((!is.null(PEM.lambda)) + (!is.null(PEM.s)) == 1) stop("Both 'PEM.s' and 'PEM.lambda' must be specified as numeric vectors or as NULL")
        if(!is.null(PEM.lambda) & !is.null(PEM.s))
        {
            if(length(PEM.lambda) != length(PEM.s)) stop("Length of starting value for 'PEM.lambda' and 'PEM.s' should be equal")
        }
        if(is.null(PEM.lambda) & is.null(PEM.s))
        {
            uniqTime <- sort(unique(Y[Y[,2]==1,1]))
            Tmax <- max(uniqTime)
            if(length(uniqTime)>20)
            {
                PEM.s <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
            }else
            {
                PEM.s <- uniqTime
            }
            PEM.lambda <- runif(length(PEM.s), -3, -2)
        }
        ##
        if(!is.null(PEM.mu_lam)){
            if(length(PEM.mu_lam) != 1) stop("Length of starting value for PEM.mu_lam should be 1")
        }
        if(is.null(PEM.mu_lam)) PEM.mu_lam <- mean(PEM.lambda)
        ##
        if(!is.null(PEM.sigSq_lam)){
            if(length(PEM.sigSq_lam) != 1) stop("Length of starting value for PEM.sigSq_lam should be 1")
        }
        if(is.null(PEM.sigSq_lam)) PEM.sigSq_lam <- ifelse(length(PEM.lambda)==1, 0.1, var(PEM.lambda))

        
        ##
        start.common <- list(beta=beta)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        start.PEM   <- list(PEM.lambda=PEM.lambda, PEM.s=PEM.s, PEM.mu_lam=PEM.mu_lam, PEM.sigSq_lam=PEM.sigSq_lam)
        
        ##
        if(!is.null(Normal.zeta)){
            if(length(Normal.zeta) != 1) stop("Length of starting value for Normal.zeta should be 1")
        }
        if(is.null(Normal.zeta)) Normal.zeta <- 1
        
        ##
        if(!is.null(DPM.class)){
            if(length(DPM.class) != J) stop(paste("Length of starting value for DPM.class must be", J))
        }
        if(is.null(DPM.class)) DPM.class <- sample(1:3, size=J, replace=TRUE)
        ##
        if(!is.null(DPM.tau)){
            if(length(DPM.tau) != 1) stop("Length of starting value for DPM.tau must be 1")
        }
        if(is.null(DPM.tau)) DPM.tau <- 0.5
        
        ##
        start.common <- list(beta=beta, V.j=V.j)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        start.PEM   <- list(PEM.lambda=PEM.lambda, PEM.s=PEM.s, PEM.mu_lam=PEM.mu_lam, PEM.sigSq_lam=PEM.sigSq_lam, K=length(PEM.s)-1)
        start.Normal    <- list(Normal.zeta=Normal.zeta)
        start.DPM    <- list(DPM.class=DPM.class, DPM.tau=DPM.tau)
        
        ##
        if(model[1] == "Weibull"){
            if(model[2] == "Normal") value <- list(common=start.common, WB=start.WB, Normal=start.Normal)
            if(model[2] == "DPM") value <- list(common=start.common, WB=start.WB, DPM=start.DPM)
        }
        if(model[1] == "PEM"){
            if(model[2] == "Normal") value <- list(common=start.common, PEM=start.PEM, Normal=start.Normal)
            if(model[2] == "DPM") value <- list(common=start.common, PEM=start.PEM, DPM=start.DPM)
        }
    }
    
    
    ## BayesID
    
    if(length(model)==2 & class(lin.pred)=="list")
    {
        print(paste("Start values are initiated for semi-competing risks ", model[2], " model...", sep = ""), cat("\n"))
        ##
        if(!is.null(cluster)) print(paste("Warning: 'cluster' is not required for ", model[2], " model so it is ignored", sep = ""))
        
        ##     
        n <- nrow(Y)
        p1 <- ncol(model.matrix(lin.pred[[1]], data=data)) - 1
        p2 <- ncol(model.matrix(lin.pred[[2]], data=data)) - 1
        p3 <- ncol(model.matrix(lin.pred[[3]], data=data)) - 1
        
        ##
        if(!is.null(beta)) stop(paste("'beta' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(PEM.lambda)) stop(paste("'PEM.lambda' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(PEM.s)) stop(paste("'PEM.s' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        ##
        if(!is.null(V.j1)) stop(paste("'V.j1' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j2)) stop(paste("'V.j2' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j3)) stop(paste("'V.j3' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        ##
        if(!is.null(MVN.SigmaV)) stop(paste("'MVN.SigmaV' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for semi-competing risks ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.class)) stop(paste("'DPM.class' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        ##
        if(!is.null(DPM.tau)) stop(paste("'DPM.tau' is for models for cluster-correlated data so it must be specified as NULL for ", model, " model.", sep = ""))
        
        ##
        if(!is.null(beta1)){
            if(length(beta1) != p1) stop(paste("Length of starting value for beta1 must be", p1))
        }
        if(is.null(beta1)) beta1 <- runif(p1, -0.1, 0.1)
        ##
        if(!is.null(beta2)){
            if(length(beta2) != p2) stop(paste("Length of starting value for beta2 must be", p2))
        }
        if(is.null(beta2)) beta2 <- runif(p2, -0.1, 0.1)
        ##
        if(!is.null(beta3)){
            if(length(beta3) != p3) stop(paste("Length of starting value for beta3 must be", p3))
        }
        if(is.null(beta3)) beta3 <- runif(p3, -0.1, 0.1)
        
        ##
        if(!is.null(theta)){
            if(length(theta) != 1) stop("Length of starting value for theta must be 1")
        }
        if(is.null(theta)) theta <- runif(1, 0.1, 1.1)
        
        if(!is.null(gamma.ji)){
            if(length(gamma.ji) != n) stop("Length of starting value for gamma.ji must be n")
        }
        if(is.null(gamma.ji)) gamma.ji <- rgamma(n, 1/theta, 1/theta)
        
        ##
        if(!is.null(WB.alpha)){
            if(length(WB.alpha) != 3) stop("Length of starting value for WB.alpha must be 3")
        }
        if(is.null(WB.alpha)) WB.alpha <- c(1, 1, 1)
        ##
        if(!is.null(WB.kappa)){
            if(length(WB.kappa) != 3) stop("Length of starting value for WB.kappa should be 3")
        }
        if(is.null(WB.kappa)) WB.kappa <- c(0.01, 0.01, 0.01)
       
        ##
        if(!is.null(PEM.s1))
        {
        	if(PEM.s1[1]== 0) stop("'min(PEM.s1)' must be greater than 0")
        }
        if((!is.null(PEM.lambda1)) + (!is.null(PEM.s1)) == 1) stop("Both 'PEM.s1' and 'PEM.lambda1' must be specified as numeric vectors or as NULL")
        if(!is.null(PEM.lambda1) & !is.null(PEM.s1))
        {
            if(length(PEM.lambda1) != length(PEM.s1)) stop("Length of starting value for 'PEM.lambda1' and 'PEM.s1' should be equal")
        }
        if(is.null(PEM.lambda1) & is.null(PEM.s1))
        {
            uniqTime <- sort(unique(Y[Y[,2]==1,1]))
            Tmax <- max(uniqTime)
            if(length(uniqTime)>20)
            {
                PEM.s1 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
            }else
            {
                PEM.s1 <- uniqTime
            }
            PEM.lambda1 <- runif(length(PEM.s1), -3, -2)
        }
        ##
        if(!is.null(PEM.s2))
        {
        	if(PEM.s2[1]== 0) stop("'min(PEM.s2)' must be greater than 0")
        }
        if((!is.null(PEM.lambda2)) + (!is.null(PEM.s2)) == 1) stop("Both 'PEM.s2' and 'PEM.lambda2' must be specified as numeric vectors or as NULL")
        if(!is.null(PEM.lambda2) & !is.null(PEM.s2))
        {
            if(length(PEM.lambda2) != length(PEM.s2)) stop("Length of starting value for 'PEM.lambda2' and 'PEM.s2' should be equal")
        }
        if(is.null(PEM.lambda2) & is.null(PEM.s2))
        {
            uniqTime <- sort(unique(Y[(Y[,2]==0) & (Y[,4]==1),3]))
            Tmax <- max(uniqTime)
            if(length(uniqTime)>20)
            {
                PEM.s2 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
            }else
            {
                PEM.s2 <- uniqTime
            }
            PEM.lambda2 <- runif(length(PEM.s2), -3, -2)
        }
        ##
        if(!is.null(PEM.s3))
        {
        	if(PEM.s3[1]== 0) stop("'min(PEM.s3)' must be greater than 0")
        }
        if((!is.null(PEM.lambda3)) + (!is.null(PEM.s3)) == 1) stop("Both 'PEM.s3' and 'PEM.lambda3' must be specified as numeric vectors or as NULL")
        if(!is.null(PEM.lambda3) & !is.null(PEM.s3))
        {
            if(length(PEM.lambda3) != length(PEM.s3)) stop("Length of starting value for 'PEM.lambda3' and 'PEM.s3' should be equal")
        }
        if(is.null(PEM.lambda3) & is.null(PEM.s3))
        {
            uniqTime <- sort(unique(Y[(Y[,2]==1) & (Y[,4]==1),3]))
            Tmax <- max(uniqTime)
            if(length(uniqTime)>20)
            {
                PEM.s3 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
            }else
            {
                PEM.s3 <- uniqTime
            }
            PEM.lambda3 <- runif(length(PEM.s3), -3, -2)
        }
        
        
        ##
        if(!is.null(PEM.mu_lam)){
            if(length(PEM.mu_lam) != 3) stop("Length of starting value for PEM.mu_lam should be 3")
        }
        if(is.null(PEM.mu_lam)) PEM.mu_lam <- c(mean(PEM.lambda1), mean(PEM.lambda2), mean(PEM.lambda3))
        ##
        if(!is.null(PEM.sigSq_lam)){
            if(length(PEM.sigSq_lam) != 3) stop("Length of starting value for PEM.sigSq_lam should be 3")
        }
        if(is.null(PEM.sigSq_lam)) PEM.sigSq_lam <- c(ifelse(length(PEM.lambda1)==1, 0.1, var(PEM.lambda1)), ifelse(length(PEM.lambda2)==1, 0.1, var(PEM.lambda2)), ifelse(length(PEM.lambda3)==1, 0.1, var(PEM.lambda3)))
        
        
        
        
        
        ##
        start.common <- list(beta1=beta1, beta2=beta2, beta3=beta3, gamma.ji=gamma.ji, theta=theta)
        start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
        start.PEM   <- list(PEM.lambda1=PEM.lambda1, PEM.lambda2=PEM.lambda2, PEM.lambda3=PEM.lambda3, PEM.s1=PEM.s1, PEM.s2=PEM.s2, PEM.s3=PEM.s3, PEM.mu_lam=PEM.mu_lam, PEM.sigSq_lam=PEM.sigSq_lam, K1=length(PEM.s1)-1, K2=length(PEM.s2)-1, K3=length(PEM.s3)-1)
        
        ##
        if(model[2] == "Weibull"){
            value <- list(common=start.common, WB=start.WB)
        }
        if(model[2] == "PEM"){
            value <- list(common=start.common, PEM=start.PEM)
        }

    }
    
   
    
    
    ## BayesIDcor
    
  	if(length(model)==3)
  	{
      print(paste("Start values are initiated for semi-competing risks ", model[2],"-", model[3], " model...", sep = ""), cat("\n"))
      ##
      if(is.null(cluster)) stop(paste("'cluster' must be given for ", model[2],"-", model[3], " model.", sep = ""))
      ##
      n <- nrow(Y)
      J <- length(unique(cluster))
      p1 <- ncol(model.matrix(lin.pred[[1]], data=data)) - 1
      p2 <- ncol(model.matrix(lin.pred[[2]], data=data)) - 1
      p3 <- ncol(model.matrix(lin.pred[[3]], data=data)) - 1

        ##
        if(!is.null(beta)) stop(paste("'beta' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))
        if(!is.null(PEM.lambda)) stop(paste("'PEM.lambda' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
        if(!is.null(PEM.s)) stop(paste("'PEM.s' is for univariate models so it must be specified as NULL for semi-competing risks ", model[2], " model.", sep = ""))
    
        if(!is.null(V.j)) stop(paste("'V.j' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))
        if(!is.null(Normal.zeta)) stop(paste("'Normal.zeta' is for univariate models so it must be specified as NULL for ", model[2],"-", model[3], " model.", sep = ""))

      ##
      if(!is.null(beta1)){
          if(length(beta1) != p1) stop(paste("Length of starting value for beta1 must be", p1))
      }
      if(is.null(beta1)) beta1 <- runif(p1, -0.1, 0.1)
      ##
      if(!is.null(beta2)){
          if(length(beta2) != p2) stop(paste("Length of starting value for beta2 must be", p2))
      }
      if(is.null(beta2)) beta2 <- runif(p2, -0.1, 0.1)
      ##
      if(!is.null(beta3)){
          if(length(beta3) != p3) stop(paste("Length of starting value for beta3 must be", p3))
      }
      if(is.null(beta3)) beta3 <- runif(p3, -0.1, 0.1)
      
      ##
      if(!is.null(V.j1)){
          if(length(V.j1) != J) stop(paste("Length of starting values for V.j1 must be", J))
      }
      if(is.null(V.j1)) V.j1 <- runif(J, -0.1, 0.1)
      ##
      if(!is.null(V.j2)){
          if(length(V.j2) != J) stop(paste("Length of starting values for V.j2 must be", J))
      }
      if(is.null(V.j2)) V.j2 <- runif(J, -0.1, 0.1)
      ##
      if(!is.null(V.j3)){
          if(length(V.j3) != J) stop(paste("Length of starting values for V.j3 must be", J))
      }
      if(is.null(V.j3)) V.j3 <- runif(J, -0.1, 0.1)
      
      ##
      if(!is.null(theta)){
          if(length(theta) != 1) stop("Length of starting value for theta must be 1")
      }
      if(is.null(theta)) theta <- runif(1, 0.1, 1.1)
      
      if(!is.null(gamma.ji)){
          if(length(gamma.ji) != n) stop("Length of starting value for gamma.ji must be n")
      }
      if(is.null(gamma.ji)) gamma.ji <- rgamma(n, 1/theta, 1/theta)
      
      ##
      if(!is.null(WB.alpha)){
          if(length(WB.alpha) != 3) stop("Length of starting value for WB.alpha must be 3")
      }
      if(is.null(WB.alpha)) WB.alpha <- c(1, 1, 1)
      ##
      if(!is.null(WB.kappa)){
          if(length(WB.kappa) != 3) stop("Length of starting value for WB.kappa should be 3")
      }
      if(is.null(WB.kappa)) WB.kappa <- c(0.01, 0.01, 0.01)
      
      ##
        if(!is.null(PEM.s1))
        {
        	if(PEM.s1[1]== 0) stop("'min(PEM.s1)' must be greater than 0")
        }
      if((!is.null(PEM.lambda1)) + (!is.null(PEM.s1)) == 1) stop("Both 'PEM.s1' and 'PEM.lambda1' must be specified as numeric vectors or as NULL")
      if(!is.null(PEM.lambda1) & !is.null(PEM.s1))
      {
          if(length(PEM.lambda1) != length(PEM.s1)) stop("Length of starting value for 'PEM.lambda1' and 'PEM.s1' should be equal")
      }
      if(is.null(PEM.lambda1) & is.null(PEM.s1))
      {
          uniqTime <- sort(unique(Y[Y[,2]==1,1]))
          Tmax <- max(uniqTime)+1
          if(length(uniqTime)>20)
          {
              PEM.s1 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
          }else
          {
              PEM.s1 <- uniqTime
          }
          PEM.lambda1 <- runif(length(PEM.s1), -3, -2)
      }
      ##
        if(!is.null(PEM.s2))
        {
        	if(PEM.s2[1]== 0) stop("'min(PEM.s2)' must be greater than 0")
        }
      if((!is.null(PEM.lambda2)) + (!is.null(PEM.s2)) == 1) stop("Both 'PEM.s2' and 'PEM.lambda2' must be specified as numeric vectors or as NULL")
      if(!is.null(PEM.lambda2) & !is.null(PEM.s2))
      {
          if(length(PEM.lambda2) != length(PEM.s2)) stop("Length of starting value for 'PEM.lambda2' and 'PEM.s2' should be equal")
      }
      if(is.null(PEM.lambda2) & is.null(PEM.s2))
      {
          uniqTime <- sort(unique(Y[(Y[,2]==0) & (Y[,4]==1),3]))
          Tmax <- max(uniqTime)+1
          if(length(uniqTime)>20)
          {
              PEM.s2 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
          }else
          {
              PEM.s2 <- uniqTime
          }
          PEM.lambda2 <- runif(length(PEM.s2), -3, -2)
      }
      ##
        if(!is.null(PEM.s3))
        {
        	if(PEM.s3[1]== 0) stop("'min(PEM.s3)' must be greater than 0")
        }
      if((!is.null(PEM.lambda3)) + (!is.null(PEM.s3)) == 1) stop("Both 'PEM.s3' and 'PEM.lambda3' must be specified as numeric vectors or as NULL")
      if(!is.null(PEM.lambda3) & !is.null(PEM.s3))
      {
          if(length(PEM.lambda3) != length(PEM.s3)) stop("Length of starting value for 'PEM.lambda3' and 'PEM.s3' should be equal")
      }
      if(is.null(PEM.lambda3) & is.null(PEM.s3))
      {
          uniqTime <- sort(unique(Y[(Y[,2]==1) & (Y[,4]==1),3]))
          Tmax <- max(uniqTime)+1
          if(length(uniqTime)>20)
          {
              PEM.s3 <- sort(c(sample(uniqTime[-length(uniqTime)], 20, replace=FALSE), Tmax))
          }else
          {
              PEM.s3 <- uniqTime
          }
          PEM.lambda3 <- runif(length(PEM.s3), -3, -2)
      }
      
      
      ##
      if(!is.null(PEM.mu_lam)){
          if(length(PEM.mu_lam) != 3) stop("Length of starting value for PEM.mu_lam should be 3")
      }
      if(is.null(PEM.mu_lam)) PEM.mu_lam <- c(mean(PEM.lambda1), mean(PEM.lambda2), mean(PEM.lambda3))
      ##
      if(!is.null(PEM.sigSq_lam)){
          if(length(PEM.sigSq_lam) != 3) stop("Length of starting value for PEM.sigSq_lam should be 3")
      }
      if(is.null(PEM.sigSq_lam)) PEM.sigSq_lam <- c(ifelse(length(PEM.lambda1)==1, 0.1, var(PEM.lambda1)), ifelse(length(PEM.lambda2)==1, 0.1, var(PEM.lambda2)), ifelse(length(PEM.lambda3)==1, 0.1, var(PEM.lambda3)))
      
      
      
      ##
      if(!is.null(MVN.SigmaV))
      {
          if(is.matrix(MVN.SigmaV) == FALSE) stop("Starting value for MVN.SigmaV must be a 3x3 matrix")
          if(is.matrix(MVN.SigmaV) == TRUE){
              if(nrow(MVN.SigmaV) != 3 | ncol(MVN.SigmaV) != 3) stop("Starting value for MVN.SigmaV must be a 3x3 matrix")
          }
      }
      if(is.null(MVN.SigmaV)) MVN.SigmaV <- diag(0.1, 3)
      
      ##
      if(!is.null(DPM.class)){
          if(length(DPM.class) != J) stop(paste("Length of starting value for DPM.class must be", J))
      }
      if(is.null(DPM.class)) DPM.class <- sample(1:3, size=J, replace=TRUE)
      ##
      if(!is.null(DPM.tau)){
          if(length(DPM.tau) != 1) stop("Length of starting value for DPM.tau must be 1")
      }
      if(is.null(DPM.tau)) DPM.tau <- 0.5
      
      ##
      start.common <- list(beta1=beta1, beta2=beta2, beta3=beta3, gamma.ji=gamma.ji, V.j1=V.j1, V.j2=V.j2, V.j3=V.j3, theta=theta)
      start.WB     <- list(WB.alpha=WB.alpha, WB.kappa=WB.kappa)
      start.PEM   <- list(PEM.lambda1=PEM.lambda1, PEM.lambda2=PEM.lambda2, PEM.lambda3=PEM.lambda3, PEM.s1=PEM.s1, PEM.s2=PEM.s2, PEM.s3=PEM.s3, PEM.mu_lam=PEM.mu_lam, PEM.sigSq_lam=PEM.sigSq_lam, K1=length(PEM.s1)-1, K2=length(PEM.s2)-1, K3=length(PEM.s3)-1)
      start.MVN    <- list(MVN.SigmaV=MVN.SigmaV)
      start.DPM    <- list(DPM.class=DPM.class, DPM.tau=DPM.tau)
      
      ##
      if(model[2] == "Weibull"){
          if(model[3] == "MVN") value <- list(common=start.common, WB=start.WB, MVN=start.MVN)
          if(model[3] == "DPM") value <- list(common=start.common, WB=start.WB, DPM=start.DPM)
      }
      if(model[2] == "PEM"){
          if(model[3] == "MVN") value <- list(common=start.common, PEM=start.PEM, MVN=start.MVN)
          if(model[3] == "DPM") value <- list(common=start.common, PEM=start.PEM, DPM=start.DPM)
      }
  }
  
  ##
  return(value)
}






