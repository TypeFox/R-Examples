


BayesSurv <- function(Y,
						lin.pred,
						data,
						 cluster=NULL,
                        model = "Weibull",
						hyperParams,
						startValues,								
						mcmc,
                        path = NULL)
{
    hyperP      <- hyperParams
    mcmcList    <- mcmc
    
    if((mcmcList$run$numReps / mcmcList$run$thin * mcmcList$run$burninPerc) %% 1 == 0)
    {
        
        # for independent univariate time-to-event data
        if(is.null(cluster))
        {
            nChain <- length(startValues)
            
            hz.type 	<- model
            
            Xmat <- model.frame(lin.pred, data=data)
            
            p <- ncol(Xmat)
            
            
            ##
            
            if(p == 0){
                survData <- Y
            }
            
            if(p > 0){
                survData <- cbind(Y, Xmat)
            }
            
            n	<- dim(survData)[1]
            
            if(!is.null(path)){
                dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
            }
            
            
            ### setting hyperparameters
            
            if(hz.type == "Weibull")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab, hyperP$WB$WB.cd))
            }
            
            if(hz.type == "PEM")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab, hyperP$PEM$PEM.alpha, 1))
            }
            
            ### mcmc setting
            
            if(hz.type == "Weibull")
            {
                mcmc <- as.vector(c(mhProp_alpha_var=mcmcList$tuning$mhProp_alpha_var, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc))
            }
            if(hz.type == "PEM")
            {
                mcmcList$tuning$s_max <- 1
                mcmc <- as.vector(c(C=mcmcList$tuning$C, delPert=mcmcList$tuning$delPert, rj.scheme = mcmcList$tuning$rj.scheme, K_max=mcmcList$tuning$K_max, s_max=mcmcList$tuning$s_max, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc))
            }
            
            
            
            
            chain = 1
            ret <- list()
            
            while(chain <= nChain){
                
                cat("chain: ", chain, "\n")
                
                nam = paste("chain", chain, sep="")
                
                temp <- startValues[[chain]]
                
                
                ### setting starting values
                
                if(hz.type == "Weibull")
                {
                    startV <- as.vector(c(beta=temp$common$beta, temp$WB$WB.alpha, temp$WB$WB.kappa))
                }
                
                if(hz.type == "PEM")
                {
                    startV <- as.vector(c(beta=temp$common$beta))
                }
                
                
                
                
                # hz.type = "Weibull"
                
                if(hz.type == "Weibull"){
                    
                    numReps     <- mcmc[2]
                    thin        <- mcmc[3]
                    burninPerc  <- mcmc[4]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    mcmcParams <- mcmc[1]
                    
                    
                    mcmcRet <- .C("BweibSurvmcmc",
                    survData 		= as.double(as.matrix(survData)),
                    n				= as.integer(n),
                    p				= as.integer(p),
                    hyperParams 	= as.double(hyperParams),
                    mcmcParams		= as.double(mcmcParams),
                    startValues 	= as.double(startV),
                    numReps			= as.integer(numReps),
                    thin			= as.integer(thin),
                    burninPerc      = as.double(burninPerc),
                    samples_beta 	= as.double(rep(0, nStore*p)),
                    samples_alpha 	= as.double(rep(0, nStore*1)),
                    samples_kappa 	= as.double(rep(0, nStore*1)),
                    samples_misc	= as.double(rep(0, p + 1)),
                    moveVec         = as.double(rep(0, numReps)))
                    
                    if(p > 0){
                        beta.p 			<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                    }
                    if(p == 0){
                        beta.p 			<- NULL
                    }
                    alpha.p 		<- matrix(mcmcRet$samples_alpha, nrow = nStore, byrow = TRUE)
                    kappa.p 		<- matrix(mcmcRet$samples_kappa, nrow = nStore, byrow = TRUE)
                    if(p > 0){
                        accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                    }
                    if(p == 0){
                        accept.beta 	<- NULL
                    }
                    accept.alpha	<- as.vector(mcmcRet$samples_misc[(p+1)])/sum(as.vector(mcmcRet$moveVec)==3)
                    if(p > 0){
                        covNames = colnames(Xmat)
                    }
                    if(p == 0){
                        covNames = NULL
                    }
                    
                    
                    ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, accept.beta = accept.beta,
                    accept.alpha = accept.alpha, covNames = covNames, hz.type = hz.type)
                    
                    
                }
                
                
                
                # hz.type = "PEM"
                
                if(hz.type == "PEM"){
                    
                    C           <- mcmc[1]
                    delPert     <- mcmc[2]
                    rj.scheme   <- mcmc[3]
                    J_max       <- mcmc[4]
                    s_max       <- max(temp$PEM$PEM.s)
                    
                    numReps     <- mcmc[6]
                    thin        <- mcmc[7]
                    burninPerc  <- mcmc[8]
                    
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    ## recommended when the unit of time is day
                    if(rj.scheme == 1){
                        s_propBI <- seq(1, s_max, 1)
                        s_propBI <- s_propBI[s_propBI < s_max]
                    }
                    ## uniquely ordered failure time
                    if(rj.scheme == 2){
                        s_propBI <- sort(unique(survData[survData[,2]==1,1]))
                        s_propBI <- s_propBI[s_propBI < s_max]
                    }
                    
                    num_s_propBI=length(s_propBI)
                    
                    time_lambda <- mcmcList$tuning$time_lambda
                    
                    nTime_lambda <- length(time_lambda)
                    
                    mcmcParams <- c(C, delPert, num_s_propBI, J_max, s_max, nTime_lambda, s_propBI, time_lambda)
                    
                    s  		<- temp$PEM$PEM.s
                    lambda  <- temp$PEM$PEM.lambda
                    
                    J=temp$PEM$K
                    mu_lam=temp$PEM$PEM.mu_lam
                    sigSq_lam=temp$PEM$PEM.sigSq_lam
                    
                    
                    J_ <- J
                    
                    startV <- as.vector(c(startV,
                    J,
                    mu_lam,
                    sigSq_lam,
                    lambda, s))
                    
                    
                    mcmcRet <- .C("BpeSurvmcmc",
                    survData            = as.double(as.matrix(survData)),
                    n                   = as.integer(n),
                    p                   = as.integer(p),
                    hyperParams         = as.double(hyperParams),
                    startValues         = as.double(startV),
                    mcmcParams          = as.double(mcmcParams),
                    numReps             = as.integer(numReps),
                    thin                = as.integer(thin),
                    burninPerc          = as.double(burninPerc),
                    samples_beta        = as.double(rep(0, nStore*p)),
                    samples_mu_lam      = as.double(rep(0, nStore*1)),
                    samples_sigSq_lam	= as.double(rep(0, nStore*1)),
                    samples_J           = as.double(rep(0, nStore*1)),
                    samples_s           = as.double(rep(0, nStore*(J_max + 1))),
                    samples_misc        = as.double(rep(0, p + 2)),
                    lambda_fin			= as.double(rep(0, nStore*nTime_lambda)),
                    moveVec             = as.double(rep(0, numReps)))
                    
                    if(p > 0){
                        beta.p 			<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                    }
                    if(p == 0){
                        beta.p 			<- NULL
                    }
                    lambda.fin 		<- matrix(mcmcRet$lambda_fin, nrow = nStore, byrow = TRUE)
                    mu_lam.p 		<- matrix(mcmcRet$samples_mu_lam, nrow = nStore, byrow = TRUE)
                    sigSq_lam.p 	<- matrix(mcmcRet$samples_sigSq_lam, nrow = nStore, byrow = TRUE)
                    K.p 			<- matrix(mcmcRet$samples_J, nrow = nStore, byrow = TRUE)
                    s.p 			<- matrix(mcmcRet$samples_s, nrow = nStore, byrow = TRUE)
                    if(p > 0){
                        accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                    }
                    if(p == 0){
                        accept.beta 	<- NULL
                    }
                    accept.BI		<- as.vector(mcmcRet$samples_misc[(p+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                    accept.DI		<- as.vector(mcmcRet$samples_misc[(p+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                    
                    if(p > 0){
                        covNames = colnames(Xmat)
                    }
                    if(p == 0){
                        covNames = NULL
                    }
                    
                    
                    ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, 
                    K.p = K.p, s.p = s.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI,
                    covNames = covNames, time_lambda = time_lambda, hz.type = hz.type)
                    
                }
                
                
                chain = chain + 1	
            }		
            
            ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, nChain = nChain)
            
            if(hz.type == "Weibull")
            {
                class(ret) <- c("Bayes", "Surv", "Ind", "WB")
            }
            if(hz.type == "PEM")
            {
                class(ret) <- c("Bayes", "Surv", "Ind", "PEM")
            }
            
            return(ret)

        }


        # for cluster-correlated univariate time-to-event data
        if(!is.null(cluster))
        {
            nChain <- length(startValues)
            
            hz.type 	<- model[1]
            re.type 	<- model[2]
            
            Xmat <- model.frame(lin.pred, data=data)
            
            p <- ncol(Xmat)
            
            
            ##
            survData <- cbind(Y, cluster, Xmat)
            
            n	<- dim(survData)[1]
            
            J	<- length(unique(cluster))
            
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(cluster == i))
            }
            
            if(!is.null(path)){
                dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
            }
            
            
            ### setting hyperparameters
            
            if(hz.type == "Weibull" & re.type == "Normal")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab, hyperP$WB$WB.cd, hyperP$Normal$Normal.ab))
            }
            if(hz.type == "Weibull" & re.type == "DPM")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab, hyperP$WB$WB.cd, 0, 1, hyperP$DPM$DPM.ab, hyperP$DPM$aTau, hyperP$DPM$bTau))
            }
            if(hz.type == "PEM" & re.type == "Normal")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab, hyperP$PEM$PEM.alpha, 1, hyperP$Normal$Normal.ab))
            }
            if(hz.type == "PEM" & re.type == "DPM")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab, hyperP$PEM$PEM.alpha, 1, 0, 1, hyperP$DPM$DPM.ab, hyperP$DPM$aTau, hyperP$DPM$bTau))
            }
            
            ### mcmc setting
            
            if(hz.type == "Weibull")
            {
                mcmc <- as.vector(c(mhProp_alpha_var=mcmcList$tuning$mhProp_alpha_var, mhProp_V_var=mcmcList$tuning$mhProp_V_var, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc, storeV=mcmcList$storage$storeV))
            }
            if(hz.type == "PEM")
            {
                mcmcList$tuning$s_max <- 1
                mcmc <- as.vector(c(C=mcmcList$tuning$C, delPert=mcmcList$tuning$delPert, rj.scheme = mcmcList$tuning$rj.scheme, K_max=mcmcList$tuning$K_max, s_max=mcmcList$tuning$s_max, mhProp_V_var=mcmcList$tuning$mhProp_V_var, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc, storeV=mcmcList$storage$storeV))
            }
            
            
            chain = 1
            ret <- list()
            
            while(chain <= nChain){
                
                cat("chain: ", chain, "\n")
                nam = paste("chain", chain, sep="")
                
                temp <- startValues[[chain]]
                
                
                ### setting starting values
                
                if(hz.type == "Weibull" & re.type == "Normal")
                {
                    startV <- as.vector(c(beta=temp$common$beta, temp$WB$WB.alpha, temp$WB$WB.kappa, V=temp$common$V.j, zeta=temp$Normal$Normal.zeta))
                }
                if(hz.type == "Weibull" & re.type == "DPM")
                {
                    startV <- as.vector(c(beta=temp$common$beta, temp$WB$WB.alpha, temp$WB$WB.kappa, V=temp$common$V.j, class=temp$DPM$DPM.class, tau=temp$DPM$DPM.tau))
                }
                if(hz.type == "PEM" & re.type == "Normal")
                {
                    startV <- as.vector(c(beta=temp$common$beta, V=temp$common$V.j, zeta=temp$Normal$Normal.zeta))
                }
                if(hz.type == "PEM" & re.type == "DPM")
                {
                    startV <- as.vector(c(beta=temp$common$beta, V=temp$common$V.j, class=temp$DPM$DPM.class, tau=temp$DPM$DPM.tau))
                }
                
                
                # hz.type = "Weibull"
                
                if(hz.type == "Weibull"){
                    
                    numReps     <- mcmc[3]
                    thin        <- mcmc[4]
                    burninPerc  <- mcmc[5]
                    storeV      <- mcmc[6]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    mcmcParams <- mcmc[1:2]
                    
                    
                    # re.type = "Normal"
                    
                    #######################################
                    ############ Weibull-Normal ###########
                    #######################################
                    
                    if(re.type == "Normal"){
                        
                        ###
                        
                        mcmcRet <- .C("BweibCorSurvmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p				= as.integer(p),
                        J				= as.integer(J),
                        nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        samples_beta 	= as.double(rep(0, nStore*p)),
                        samples_alpha 	= as.double(rep(0, nStore*1)),
                        samples_kappa 	= as.double(rep(0, nStore*1)),
                        samples_V		= as.double(rep(0, nStore*J)),
                        samples_zeta	= as.double(rep(0, nStore*1)),
                        samples_misc	= as.double(rep(0, p+1+J)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        
                        
                        if(p > 0){
                            beta.p 		<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                        }
                        if(p == 0){
                            beta.p 		<- NULL
                        }
                        
                        
                        alpha.p 		<- matrix(mcmcRet$samples_alpha, nrow = nStore, byrow = TRUE)
                        kappa.p 		<- matrix(mcmcRet$samples_kappa, nrow = nStore, byrow = TRUE)
                        V.p             <- matrix(mcmcRet$samples_V, nrow = nStore, byrow = TRUE)
                        zeta.p 			<- matrix(mcmcRet$samples_zeta, nrow = nStore, byrow = TRUE)
                        
                        if(p > 0){
                            accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                        }
                        if(p == 0){
                            accept.beta 	<- NULL
                        }
                        
                        accept.alpha	<- as.vector(mcmcRet$samples_misc[p+1])/sum(as.vector(mcmcRet$moveVec)==2)
                        accept.V 	<- as.vector(mcmcRet$samples_misc[(p+1+1):(p+1+J)])/sum(as.vector(mcmcRet$moveVec)==4)
                        
                        Vsummary <- as.matrix(apply(V.p, 2, summary))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
                        rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
                        
                        if(storeV == TRUE & !is.null(path))
                        {
                            save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
                        }
                        
                        if(p > 0){
                            covNames = colnames(Xmat)
                        }
                        if(p == 0){
                            covNames = NULL
                        }
                        
                        ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, zeta.p = zeta.p, accept.beta = accept.beta, accept.alpha = accept.alpha, accept.V = accept.V, Vsum = Vsummary, covNames = covNames)
                        
                        
                    } ## end: if Weibull-Normal
                    
                    # re.type = "DPM"
                    
                    #################################################
                    ############ Weibull-DPM (univariate) ###########
                    #################################################
                    
                    if(re.type == "DPM"){
                        
                        mcmcRet <- .C("BweibDpCorSurvmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p				= as.integer(p),
                        J				= as.integer(J),
                        nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        samples_beta 	= as.double(rep(0, nStore*p)),
                        samples_alpha 	= as.double(rep(0, nStore*1)),
                        samples_kappa 	= as.double(rep(0, nStore*1)),
                        samples_V		= as.double(rep(0, nStore*J)),
                        samples_c		= as.double(rep(0, nStore*J)),
                        samples_mu		= as.double(rep(0, nStore*J)),
                        samples_zeta	= as.double(rep(0, nStore*J)),
                        samples_tau     = as.double(rep(0, nStore*1)),
                        samples_misc	= as.double(rep(0, p+1+J)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        if(p > 0){
                            beta.p 		<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                        }
                        if(p == 0){
                            beta.p 		<- NULL
                        }
                        
                        
                        alpha.p 		<- matrix(mcmcRet$samples_alpha, nrow = nStore, byrow = TRUE)
                        kappa.p 		<- matrix(mcmcRet$samples_kappa, nrow = nStore, byrow = TRUE)
                        V.p             <- matrix(mcmcRet$samples_V, nrow = nStore, byrow = TRUE)
                        c.p				<- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                        mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                        zeta.p          <- matrix(mcmcRet$samples_zeta, nrow = nStore, byrow = TRUE)
                        tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                        
                        if(p > 0){
                            accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                        }
                        if(p == 0){
                            accept.beta 	<- NULL
                        }
                        
                        accept.alpha	<- as.vector(mcmcRet$samples_misc[p+1])/sum(as.vector(mcmcRet$moveVec)==2)
                        accept.V 	<- as.vector(mcmcRet$samples_misc[(p+1+1):(p+1+J)])/sum(as.vector(mcmcRet$moveVec)==4)
                        
                        
                        Vsummary <- as.matrix(apply(V.p, 2, summary))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
                        rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
                        
                        if(storeV == TRUE & !is.null(path))
                        {
                            save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
                        }
                        
                        if(p > 0){
                            covNames = colnames(Xmat)
                        }
                        if(p == 0){
                            covNames = NULL
                        }
                        
                        ret[[nam]] <- list(beta.p = beta.p, alpha.p = alpha.p, kappa.p = kappa.p, class.p = c.p, mu.p = mu.p, zeta.p = zeta.p, tau.p = tau.p, accept.beta = accept.beta, accept.alpha = accept.alpha, accept.V = accept.V, covNames = covNames, Vsum = Vsummary)
                        
                    } ## end: if Weibull-DPM
                    
                } ## end: if Weibull
                
                
                # hz.type = "PEM"
                
                if(hz.type == "PEM"){
                    
                    C           <- mcmc[1]
                    delPert     <- mcmc[2]
                    rj.scheme   <- mcmc[3]
                    K_max       <- mcmc[4]
                    s_max       <- max(temp$PEM$PEM.s)
                    mhProp_V_var <- mcmc[6]
                    numReps     <- mcmc[7]
                    thin        <- mcmc[8]
                    burninPerc  <- mcmc[9]
                    storeV      <- mcmc[10]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    ## recommended when the unit of time is day
                    if(rj.scheme == 1){
                        s_propBI <- seq(1, s_max, 1)
                        s_propBI <- s_propBI[s_propBI < s_max]
                    }
                    ## uniquely ordered failure time
                    if(rj.scheme == 2){
                        s_propBI <- sort(unique(survData[survData[,2]==1,1]))
                        s_propBI <- s_propBI[s_propBI < s_max]
                    }
                    
                    num_s_propBI=length(s_propBI)
                    
                    time_lambda <- mcmcList$tuning$time_lambda
                    
                    nTime_lambda <- length(time_lambda)
                    
                    mcmcParams <- c(C, delPert, num_s_propBI, K_max, s_max, nTime_lambda, s_propBI, time_lambda, mhProp_V_var)
                    
                    
                    s  		<- temp$PEM$PEM.s
                    lambda <- temp$PEM$PEM.lambda
                    
                    K=temp$PEM$K
                    mu_lam=temp$PEM$PEM.mu_lam
                    sigSq_lam=temp$PEM$PEM.sigSq_lam
                    
                    
                    
                    # re.type = "Normal"
                    
                    ###################################
                    ############ PEM-Normal ###########
                    ###################################
                    
                    if(re.type == "Normal"){
                        
                        ###
                        K_ <- K
                        
                        startV <- as.vector(c(startV[1:(p)],
                        K,
                        mu_lam,
                        sigSq_lam,
                        lambda, s,
                        startV[(p+1):length(startV)]))
                        
                        
                        mcmcRet <- .C("BpeMvnCorSurvmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p				= as.integer(p),
                        J				= as.integer(J),
                        nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        samples_beta 	= as.double(rep(0, nStore*p)),
                        samples_mu_lam  = as.double(rep(0, nStore*1)),
                        samples_sigSq_lam  = as.double(rep(0, nStore*1)),
                        samples_K          = as.double(rep(0, nStore*1)),
                        samples_s          = as.double(rep(0, nStore*(K_max + 1))),
                        samples_V		= as.double(rep(0, nStore*J)),
                        samples_zeta	= as.double(rep(0, nStore*1)),
                        samples_misc	= as.double(rep(0, p+2+J)),
                        lambda_fin     = as.double(rep(0, nStore*nTime_lambda)),
                        dev     		= as.double(rep(0, nStore*1)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        
                        
                        if(p > 0){
                            beta.p 		<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                        }
                        if(p == 0){
                            beta.p 		<- NULL
                        }
                        
                        lambda.fin 	<- matrix(mcmcRet$lambda_fin, nrow = nStore, byrow = TRUE)
                        
                        mu_lam.p 		<- matrix(mcmcRet$samples_mu_lam, nrow = nStore, byrow = TRUE)
                        sigSq_lam.p 	<- matrix(mcmcRet$samples_sigSq_lam, nrow = nStore, byrow = TRUE)
                        K.p 			<- matrix(mcmcRet$samples_K, nrow = nStore, byrow = TRUE)
                        s.p 			<- matrix(mcmcRet$samples_s, nrow = nStore, byrow = TRUE)
                        V.p             <- matrix(mcmcRet$samples_V, nrow = nStore, byrow = TRUE)
                        zeta.p 			<- matrix(mcmcRet$samples_zeta, nrow = nStore, byrow = TRUE)
                        
                        if(p > 0){
                            accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                        }
                        if(p == 0){
                            accept.beta 	<- NULL
                        }
                        
                        accept.BI		<- as.vector(mcmcRet$samples_misc[(p)+1])/sum(as.vector(mcmcRet$moveVec)==4)
                        accept.DI		<- as.vector(mcmcRet$samples_misc[(p)+2])/sum(as.vector(mcmcRet$moveVec)==5)
                        accept.V 	<- as.vector(mcmcRet$samples_misc[(p+3):(p+2+J)])/sum(as.vector(mcmcRet$moveVec)==6)
                        
                        dev.p 			<- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                        
                        Vsummary <- as.matrix(apply(V.p, 2, summary))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
                        rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
                        
                        if(storeV == TRUE & !is.null(path))
                        {
                            save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
                        }
                        
                        if(p > 0){
                            covNames = colnames(Xmat)
                        }
                        if(p == 0){
                            covNames = NULL
                        }
  
                        ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, K.p = K.p, s.p = s.p, zeta.p = zeta.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI, time_lambda = time_lambda, accept.V = accept.V, covNames = covNames, Vsum = Vsummary)
                        
                    }   ## end: if PEM-Normal
                    
                    # re.type = "DPM"
                    
                    ###################################
                    ############ PEM-DPM ###########
                    ###################################
                    
                    if(re.type == "DPM"){
                        
                        ###
                        K_ <- K
                        
                        startV <- as.vector(c(startV[1:(p)],
                        K,
                        mu_lam,
                        sigSq_lam,
                        lambda, s,
                        startV[(p+1):length(startV)]))
                        
                        
                        
                        mcmcRet <- .C("BpeDpCorSurvmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p				= as.integer(p),
                        J				= as.integer(J),
                        nj				= as.double(nj),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        samples_beta 	= as.double(rep(0, nStore*p)),
                        samples_mu_lam  = as.double(rep(0, nStore*1)),
                        samples_sigSq_lam  = as.double(rep(0, nStore*1)),
                        samples_K          = as.double(rep(0, nStore*1)),
                        samples_s          = as.double(rep(0, nStore*(K_max + 1))),
                        samples_V		= as.double(rep(0, nStore*J)),
                        samples_c		= as.double(rep(0, nStore*J)),
                        samples_mu		= as.double(rep(0, nStore*J)),
                        samples_zeta	= as.double(rep(0, nStore*J)),
                        samples_tau     = as.double(rep(0, nStore*1)),
                        samples_misc	= as.double(rep(0, p+2+J)),
                        lambda_fin     = as.double(rep(0, nStore*nTime_lambda)),
                        dev     		= as.double(rep(0, nStore*1)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        
                        
                        if(p > 0){
                            beta.p 		<- matrix(mcmcRet$samples_beta, nrow = nStore, byrow = TRUE)
                        }
                        if(p == 0){
                            beta.p 		<- NULL
                        }
                        
                        lambda.fin 	<- matrix(mcmcRet$lambda_fin, nrow = nStore, byrow = TRUE)
                        
                        mu_lam.p 		<- matrix(mcmcRet$samples_mu_lam, nrow = nStore, byrow = TRUE)
                        sigSq_lam.p 	<- matrix(mcmcRet$samples_sigSq_lam, nrow = nStore, byrow = TRUE)
                        K.p 			<- matrix(mcmcRet$samples_K, nrow = nStore, byrow = TRUE)
                        s.p 			<- matrix(mcmcRet$samples_s, nrow = nStore, byrow = TRUE)
                        V.p             <- matrix(mcmcRet$samples_V, nrow = nStore, byrow = TRUE)
                        c.p				<- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                        mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                        zeta.p          <- matrix(mcmcRet$samples_zeta, nrow = nStore, byrow = TRUE)
                        tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                        
                        if(p > 0){
                            accept.beta 	<- as.vector(mcmcRet$samples_misc[1:p])/sum(as.vector(mcmcRet$moveVec)==1)*p
                        }
                        if(p == 0){
                            accept.beta 	<- NULL
                        }
                        
                        accept.BI		<- as.vector(mcmcRet$samples_misc[(p)+1])/sum(as.vector(mcmcRet$moveVec)==4)
                        accept.DI		<- as.vector(mcmcRet$samples_misc[(p)+2])/sum(as.vector(mcmcRet$moveVec)==5)
                        accept.V 	<- as.vector(mcmcRet$samples_misc[(p+3):(p+2+J)])/sum(as.vector(mcmcRet$moveVec)==6)
                        
                        dev.p 			<- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                        
                        Vsummary <- as.matrix(apply(V.p, 2, summary))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.975))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, quantile, prob = 0.025))
                        Vsummary <- rbind(Vsummary, apply(V.p, 2, sd))
                        rownames(Vsummary)[7:9] <- c("0.975", "0.025", "sd")
                        
                        if(storeV == TRUE)
                        {
                            save(V.p, file = paste(path, "VPch", chain, ".RData", sep = ""))
                        }
                        
                        if(p > 0){
                            covNames = colnames(Xmat)
                        }
                        if(p == 0){
                            covNames = NULL
                        }
                        
                        
                        ret[[nam]] <- list(beta.p = beta.p, lambda.fin = lambda.fin, mu_lam.p = mu_lam.p, sigSq_lam.p = sigSq_lam.p, K.p = K.p, s.p = s.p, class.p = c.p, mu.p = mu.p, zeta.p = zeta.p, tau.p = tau.p, accept.beta = accept.beta, accept.BI = accept.BI, accept.DI = accept.DI, time_lambda = time_lambda, accept.V = accept.V, covNames = covNames, Vsum = Vsummary)
                        
                        
                    }   ## end: if PEM-DPM
                    
                    
                } ## end: if PEM
                
                
                
                chain = chain + 1	
                
            }## end: while(chain <= nChain)
            
            
            
            ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, re.type = re.type, nChain = nChain)
            
            if(hz.type == "Weibull")
            {
                if(re.type == "Normal")
                {
                    class(ret) <- c("Bayes", "Surv", "Cor", "WB", "Normal")
                }
                if(re.type == "DPM")
                {
                    class(ret) <- c("Bayes", "Surv", "Cor", "WB", "DPM")
                }
            }
            if(hz.type == "PEM")
            {
                if(re.type == "Normal")
                {
                    class(ret) <- c("Bayes", "Surv", "Cor", "PEM", "Normal")
                }
                if(re.type == "DPM")
                {
                    class(ret) <- c("Bayes", "Surv", "Cor", "PEM", "DPM")
                }
            }
            return(ret)
        }



    }
    else{
    	print(" (numReps * burninPerc) must be divisible by (thin)")
    }
    
	

}





























