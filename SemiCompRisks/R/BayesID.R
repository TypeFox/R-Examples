


BayesID <- function(Y,
                    lin.pred,
                    data,
                    cluster = NULL,
                    model = c("semi-Markov", "Weibull"),
                    hyperParams,
                    startValues,
                    mcmc,
                    path = NULL)
{
    hyperP      <- hyperParams
    mcmcList    <- mcmc
    
    if((mcmcList$run$numReps / mcmcList$run$thin * mcmcList$run$burninPerc) %% 1 == 0)
    {
        ## for independent semi-competing risks data
        
        if(is.null(cluster))
        {
            nChain <- length(startValues)
            
            model_h3    <- model[1]
            hz.type     <- model[2]
            
            Xmat1 <- model.frame(lin.pred[[1]], data=data)
            Xmat2 <- model.frame(lin.pred[[2]], data=data)
            Xmat3 <- model.frame(lin.pred[[3]], data=data)
            
            p1 <- ncol(Xmat1)
            p2 <- ncol(Xmat2)
            p3 <- ncol(Xmat3)
            
            nCov <- c(p1, p2, p3)
            
            ##
            
            if(p1 == 0 & p2 == 0 & p3 == 0){
                survData <- Y
            }
            if(p1 == 0 & p2 > 0 & p3 > 0){
                survData <- cbind(Y, Xmat2, Xmat3)
            }
            if(p1 > 0 & p2 == 0 & p3 > 0){
                survData <- cbind(Y, Xmat1, Xmat3)
            }
            if(p1 > 0 & p2 > 0 & p3 == 0){
                survData <- cbind(Y, Xmat1, Xmat2)
            }
            if(p1 > 0 & p2 == 0 & p3 == 0){
                survData <- cbind(Y, Xmat1)
            }
            if(p1 == 0 & p2 > 0 & p3 == 0){
                survData <- cbind(Y, Xmat2)
            }
            if(p1 == 0 & p2 == 0 & p3 > 0){
                survData <- cbind(Y, Xmat3)
            }
            if(p1 > 0 & p2 > 0 & p3 > 0){
                survData <- cbind(Y, Xmat1, Xmat2, Xmat3)
            }
            
            n	<- dim(survData)[1]
            
            #if(class(startValues) == "list" & length(startValues) == nChain){
            
            if(!is.null(path)){
                dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
            }
            
            ### setting hyperparameters
            
            if(hz.type == "Weibull")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab1, hyperP$WB$WB.ab2, hyperP$WB$WB.ab3, hyperP$WB$WB.cd1, hyperP$WB$WB.cd2, hyperP$WB$WB.cd3, hyperP$theta))
            }
            
            if(hz.type == "PEM")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab1, hyperP$PEM$PEM.ab2, hyperP$PEM$PEM.ab3, hyperP$PEM$PEM.alpha1, hyperP$PEM$PEM.alpha2, hyperP$PEM$PEM.alpha3, 1, 1, 1, hyperP$theta))
            }
            
            
            
            ### mcmc setting
            
            if(hz.type == "Weibull")
            {
                mcmc <- as.vector(c(mhProp_alpha1_var=mcmcList$tuning$mhProp_alphag_var[1], mhProp_alpha2_var=mcmcList$tuning$mhProp_alphag_var[2], mhProp_alpha3_var=mcmcList$tuning$mhProp_alphag_var[3], mhProp_theta_var=mcmcList$tuning$mhProp_theta_var, nGam_save=mcmcList$storage$nGam_save, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc))
            }
            if(hz.type == "PEM")
            {
                mcmcList$tuning$sg_max <- c(1,1,1)
                mcmc <- as.vector(c(C1=mcmcList$tuning$Cg[1], C2=mcmcList$tuning$Cg[2], C3=mcmcList$tuning$Cg[3], delPert1=mcmcList$tuning$delPertg[1], delPert2=mcmcList$tuning$delPertg[2], delPert3=mcmcList$tuning$delPertg[3], rj.scheme = mcmcList$tuning$rj.scheme, K1_max=mcmcList$tuning$Kg_max[1], K2_max=mcmcList$tuning$Kg_max[2], K3_max=mcmcList$tuning$Kg_max[3], s1_max=mcmcList$tuning$sg_max[1], s2_max=mcmcList$tuning$sg_max[2], s3_max=mcmcList$tuning$sg_max[3], mhProp_theta_var=mcmcList$tuning$mhProp_theta_var, nGam_save=mcmcList$storage$nGam_save, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc))
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
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, temp$WB$WB.alpha, temp$WB$WB.kappa, theta=temp$common$theta, gamma=temp$common$gamma.ji))
                }
                if(hz.type == "PEM")
                {
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, theta=temp$common$theta, gamma=temp$common$gamma.ji))
                }
                
                
                # hz.type = "Weibull"
                
                if(hz.type == "Weibull"){
                    
                    nGam_save   <- mcmc[5]
                    numReps     <- mcmc[6]
                    thin        <- mcmc[7]
                    burninPerc  <- mcmc[8]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    mcmcParams <- mcmc[1:4]
                    
                    # model_h3 = "Markov"
                    
                    if(model_h3 == "Markov"){
                        
                        ###
                        
                        mcmcRet <- .C("BweibScrmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p1				= as.integer(p1),
                        p2				= as.integer(p2),
                        p3				= as.integer(p3),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        nGam_save		= as.integer(nGam_save),
                        samples_beta1 	= as.double(rep(0, nStore*p1)),
                        samples_beta2 	= as.double(rep(0, nStore*p2)),
                        samples_beta3 	= as.double(rep(0, nStore*p3)),
                        samples_alpha1 	= as.double(rep(0, nStore*1)),
                        samples_alpha2 	= as.double(rep(0, nStore*1)),
                        samples_alpha3 	= as.double(rep(0, nStore*1)),
                        samples_kappa1 	= as.double(rep(0, nStore*1)),
                        samples_kappa2 	= as.double(rep(0, nStore*1)),
                        samples_kappa3 	= as.double(rep(0, nStore*1)),
                        samples_theta 	= as.double(rep(0, nStore*1)),
                        samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                        samples_misc	= as.double(rep(0, (p1+p2+p3+4))),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        if(p1 > 0){
                            beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                        }
                        if(p1 == 0){
                            beta1.p 		<- NULL
                        }
                        if(p2 > 0){
                            beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                        }
                        if(p2 == 0){
                            beta2.p 		<- NULL
                        }
                        if(p3 > 0){
                            beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                        }
                        if(p3 == 0){
                            beta3.p 		<- NULL
                        }
                        
                        alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                        alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                        alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                        kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                        kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                        kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                        theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                        gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                        
                        if(p1 > 0){
                            accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                        }
                        if(p1 == 0){
                            accept.beta1 	<- NULL
                        }
                        if(p2 > 0){
                            accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                        }
                        if(p2 == 0){
                            accept.beta2 	<- NULL
                        }
                        if(p3 > 0){
                            accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                        }
                        if(p3 == 0){
                            accept.beta3 	<- NULL
                        }
                        
                        accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                        accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                        accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                        accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                        
                        if(nGam_save > 0 & !is.null(path)){
                            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                        }
                        
                        if(p1 > 0){
                            covNames1 = colnames(Xmat1)
                        }
                        if(p1 == 0){
                            covNames1 = NULL
                        }
                        if(p2 > 0){
                            covNames2 = colnames(Xmat2)
                        }
                        if(p2 == 0){
                            covNames2 = NULL
                        }
                        if(p3 > 0){
                            covNames3 = colnames(Xmat3)
                        }
                        if(p3 == 0){
                            covNames3 = NULL
                        }
                        
                        
                        ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3)
                    } # if(model_h3 == "Markov")
                    
                    # model_h3 = "semi-Markov"
                    
                    if(model_h3 == "semi-Markov"){
            
                        ###
                        
                        mcmcRet <- .C("BweibScrSMmcmc",
                        survData 		= as.double(as.matrix(survData)),
                        n				= as.integer(n),
                        p1				= as.integer(p1),
                        p2				= as.integer(p2),
                        p3				= as.integer(p3),
                        hyperParams 	= as.double(hyperParams),
                        mcmcParams		= as.double(mcmcParams),
                        startValues 	= as.double(startV),
                        numReps			= as.integer(numReps),
                        thin			= as.integer(thin),
                        burninPerc      = as.double(burninPerc),
                        nGam_save		= as.integer(nGam_save),
                        samples_beta1 	= as.double(rep(0, nStore*p1)),
                        samples_beta2 	= as.double(rep(0, nStore*p2)),
                        samples_beta3 	= as.double(rep(0, nStore*p3)),
                        samples_alpha1 	= as.double(rep(0, nStore*1)),
                        samples_alpha2 	= as.double(rep(0, nStore*1)),
                        samples_alpha3 	= as.double(rep(0, nStore*1)),
                        samples_kappa1 	= as.double(rep(0, nStore*1)),
                        samples_kappa2 	= as.double(rep(0, nStore*1)),
                        samples_kappa3 	= as.double(rep(0, nStore*1)),
                        samples_theta 	= as.double(rep(0, nStore*1)),
                        samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                        samples_misc	= as.double(rep(0, (p1+p2+p3+4))),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        if(p1 > 0){
                            beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                        }
                        if(p1 == 0){
                            beta1.p 		<- NULL
                        }
                        if(p2 > 0){
                            beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                        }
                        if(p2 == 0){
                            beta2.p 		<- NULL
                        }
                        if(p3 > 0){
                            beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                        }
                        if(p3 == 0){
                            beta3.p 		<- NULL
                        }
                        
                        alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                        alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                        alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                        kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                        kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                        kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                        theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                        gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                        
                        if(p1 > 0){
                            accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                        }
                        if(p1 == 0){
                            accept.beta1 	<- NULL
                        }
                        if(p2 > 0){
                            accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                        }
                        if(p2 == 0){
                            accept.beta2 	<- NULL
                        }
                        if(p3 > 0){
                            accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                        }
                        if(p3 == 0){
                            accept.beta3 	<- NULL
                        }
                        
                        accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                        accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                        accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                        accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                        
                        if(nGam_save > 0 & !is.null(path)){
                            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                        }
                        
                        if(p1 > 0){
                            covNames1 = colnames(Xmat1)
                        }
                        if(p1 == 0){
                            covNames1 = NULL
                        }
                        if(p2 > 0){
                            covNames2 = colnames(Xmat2)
                        }
                        if(p2 == 0){
                            covNames2 = NULL
                        }
                        if(p3 > 0){
                            covNames3 = colnames(Xmat3)
                        }
                        if(p3 == 0){
                            covNames3 = NULL
                        }
                        
                        ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3)
                    } # if(model_h3 == "semi-Markov")
                    
                    
                } # if(hz.type == "Weibull")
                
                # hz.type = "PEM"
                
                if(hz.type == "PEM"){
                    
                    C1 <- mcmc[1]
                    C2 <- mcmc[2]
                    C3 <- mcmc[3]
                    delPert1 <- mcmc[4]
                    delPert2 <- mcmc[5]
                    delPert3 <- mcmc[6]
                    rj.scheme <- mcmc[7]
                    K1_max <- mcmc[8]
                    K2_max <- mcmc[9]
                    K3_max <- mcmc[10]
                    s1_max <- max(temp$PEM$PEM.s1)
                    s2_max <- max(temp$PEM$PEM.s2)
                    s3_max <- max(temp$PEM$PEM.s3)
                    mhProp_theta_var <- mcmc[14]
                    nGam_save   <- mcmc[15]
                    numReps     <- mcmc[16]
                    thin        <- mcmc[17]
                    burninPerc  <- mcmc[18]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    ## recommended when the unit of time is day
                    if(rj.scheme == 1){
                        s_propBI1 <- seq(1, s1_max, 1)
                        s_propBI1 <- s_propBI1[s_propBI1 < s1_max]
                        s_propBI2 <- seq(1, s2_max, 1)
                        s_propBI2 <- s_propBI2[s_propBI2 < s2_max]
                        s_propBI3 <- seq(1, s3_max, 1)
                        s_propBI3 <- s_propBI3[s_propBI3 < s3_max]
                    }
                    ## uniquely ordered failure time
                    if(rj.scheme == 2){
                        s_propBI1 <- sort(unique(survData[survData[,2]==1,1]))
                        s_propBI1 <- s_propBI1[s_propBI1 < s1_max]
                        s_propBI2 <- sort(unique(survData[(survData[,2]==0) & (survData[,4]==1),3]))
                        s_propBI2 <- s_propBI2[s_propBI2 < s2_max]
                        s_propBI3 <- sort(unique(survData[(survData[,2]==1) & (survData[,4]==1),3]))
                        s_propBI3 <- s_propBI3[s_propBI3 < s3_max]
                    }
                    
                    
                    num_s_propBI1=length(s_propBI1)
                    num_s_propBI2=length(s_propBI2)
                    num_s_propBI3=length(s_propBI3)
                    
                    time_lambda1 <- mcmcList$tuning$time_lambda1
                    time_lambda2 <- mcmcList$tuning$time_lambda2
                    time_lambda3 <- mcmcList$tuning$time_lambda3
                    
                    nTime_lambda1 <- length(time_lambda1)
                    nTime_lambda2 <- length(time_lambda2)
                    nTime_lambda3 <- length(time_lambda3)
                    
                    mcmcParams <- c(C1, C2, C3, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2, num_s_propBI3, K1_max, K2_max, K3_max, s1_max, s2_max, s3_max, nTime_lambda1, nTime_lambda2, nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1, time_lambda2, time_lambda3, mhProp_theta_var)
                    
                    
                    s1  		<- temp$PEM$PEM.s1
                    s2			<- temp$PEM$PEM.s2
                    s3			<- temp$PEM$PEM.s3
                    lambda1 <- temp$PEM$PEM.lambda1
                    lambda2 <- temp$PEM$PEM.lambda2
                    lambda3 <- temp$PEM$PEM.lambda3
                    
                    K1=temp$PEM$K1
                    K2=temp$PEM$K2
                    K3=temp$PEM$K3
                    mu_lam1=temp$PEM$PEM.mu_lam[1]
                    mu_lam2=temp$PEM$PEM.mu_lam[2]
                    mu_lam3=temp$PEM$PEM.mu_lam[3]
                    sigSq_lam1=temp$PEM$PEM.sigSq_lam[1]
                    sigSq_lam2=temp$PEM$PEM.sigSq_lam[2]
                    sigSq_lam3=temp$PEM$PEM.sigSq_lam[3]
                    
                    
                    # model_h3 = "Markov"
                    
                    if(model_h3 == "Markov"){
                        
                        ###
                        
                        J_1 <- K1
                        J_2 <- K2
                        J_3 <- K3
                        
                        theta <- startV[(p1+p2+p3+1)]
                        gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                        
                        startV <- as.vector(c(startV[1:(p1+p2+p3)],
                        K1, K2, K3,
                        mu_lam1, mu_lam2, mu_lam3,
                        sigSq_lam1, sigSq_lam2, sigSq_lam3,
                        theta, gamma,
                        lambda1, lambda2, lambda3, s1, s2, s3))
                        
                        
                        mcmcRet <- .C("BpeScrmcmc",
                        survData            = as.double(as.matrix(survData)),
                        n                   = as.integer(n),
                        p1					= as.integer(p1),
                        p2					= as.integer(p2),
                        p3					= as.integer(p3),
                        hyperParams         = as.double(hyperParams),
                        startValues         = as.double(startV),
                        mcmcParams          = as.double(mcmcParams),
                        numReps             = as.integer(numReps),
                        thin                = as.integer(thin),
                        burninPerc          = as.double(burninPerc),
                        nGam_save			= as.integer(nGam_save),
                        samples_beta1 		= as.double(rep(0, nStore*p1)),
                        samples_beta2 		= as.double(rep(0, nStore*p2)),
                        samples_beta3 		= as.double(rep(0, nStore*p3)),
                        samples_mu_lam1     = as.double(rep(0, nStore*1)),
                        samples_mu_lam2     = as.double(rep(0, nStore*1)),
                        samples_mu_lam3     = as.double(rep(0, nStore*1)),
                        samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                        samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                        samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                        samples_J1          = as.double(rep(0, nStore*1)),
                        samples_J2          = as.double(rep(0, nStore*1)),
                        samples_J3          = as.double(rep(0, nStore*1)),
                        samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                        samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                        samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                        samples_theta 		= as.double(rep(0, nStore*1)),
                        samples_gamma 		= as.double(rep(0, nStore*nGam_save)),
                        samples_misc        = as.double(rep(0, p1 + p2 + p3 + 7)),
                        lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                        lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                        lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        if(p1 > 0){
                            beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                        }
                        if(p1 == 0){
                            beta1.p 		<- NULL
                        }
                        if(p2 > 0){
                            beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                        }
                        if(p2 == 0){
                            beta2.p 		<- NULL
                        }
                        if(p3 > 0){
                            beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                        }
                        if(p3 == 0){
                            beta3.p 		<- NULL
                        }
                        
                        lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                        lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                        lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                        gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                        theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                        mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                        mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                        mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                        sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                        sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                        sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                        K1.p 			<- matrix(mcmcRet$samples_J1, nrow = nStore, byrow = TRUE)
                        K2.p 			<- matrix(mcmcRet$samples_J2, nrow = nStore, byrow = TRUE)
                        K3.p 			<- matrix(mcmcRet$samples_J3, nrow = nStore, byrow = TRUE)
                        s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                        s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                        s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                        
                        if(p1 > 0){
                            accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                        }
                        if(p1 == 0){
                            accept.beta1 	<- NULL
                        }
                        if(p2 > 0){
                            accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                        }
                        if(p2 == 0){
                            accept.beta2 	<- NULL
                        }
                        if(p3 > 0){
                            accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                        }
                        if(p3 == 0){
                            accept.beta3 	<- NULL
                        }
                        accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                        accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                        accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                        accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                        accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                        accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                        accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+7])/sum(as.vector(mcmcRet$moveVec)==11)
                        
                        if(nGam_save > 0 & !is.null(path)){
                            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                        }
                        
                        if(p1 > 0){
                            covNames1 = colnames(Xmat1)
                        }
                        if(p1 == 0){
                            covNames1 = NULL
                        }
                        
                        if(p2 > 0){
                            covNames2 = colnames(Xmat2)
                        }
                        if(p2 == 0){
                            covNames2 = NULL
                        }
                        if(p3 > 0){
                            covNames3 = colnames(Xmat3)
                        }
                        if(p3 == 0){
                            covNames3 = NULL
                        }
                        

                        
                        
                        ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, theta.p = theta.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3)
                        
                    } #if(model_h3 == "Markov")
                    
                    # model_h3 = "semi-Markov"
                    
                    if(model_h3 == "semi-Markov"){
                        
                        ###
                        
                        J_1 <- K1
                        J_2 <- K2
                        J_3 <- K3
                        
                        theta <- startV[(p1+p2+p3+1)]
                        gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                        
                        startV <- as.vector(c(startV[1:(p1+p2+p3)],
                        K1, K2, K3,
                        mu_lam1, mu_lam2, mu_lam3,
                        sigSq_lam1, sigSq_lam2, sigSq_lam3,
                        theta, gamma,
                        lambda1, lambda2, lambda3, s1, s2, s3))
                        
                        
                        
                        mcmcRet <- .C("BpeScrSMmcmc",
                        survData            = as.double(as.matrix(survData)),
                        n                   = as.integer(n),
                        p1					= as.integer(p1),
                        p2					= as.integer(p2),
                        p3					= as.integer(p3),
                        hyperParams         = as.double(hyperParams),
                        startValues         = as.double(startV),
                        mcmcParams          = as.double(mcmcParams),
                        numReps             = as.integer(numReps),
                        thin                = as.integer(thin),
                        burninPerc          = as.double(burninPerc),
                        nGam_save			= as.integer(nGam_save),
                        samples_beta1 		= as.double(rep(0, nStore*p1)),
                        samples_beta2 		= as.double(rep(0, nStore*p2)),
                        samples_beta3 		= as.double(rep(0, nStore*p3)),
                        samples_mu_lam1     = as.double(rep(0, nStore*1)),
                        samples_mu_lam2     = as.double(rep(0, nStore*1)),
                        samples_mu_lam3     = as.double(rep(0, nStore*1)),
                        samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                        samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                        samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                        samples_J1          = as.double(rep(0, nStore*1)),
                        samples_J2          = as.double(rep(0, nStore*1)),
                        samples_J3          = as.double(rep(0, nStore*1)),
                        samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                        samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                        samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                        samples_theta 		= as.double(rep(0, nStore*1)),
                        samples_gamma 		= as.double(rep(0, nStore*nGam_save)),
                        samples_misc        = as.double(rep(0, p1 + p2 + p3 + 7)),
                        lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                        lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                        lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                        moveVec             = as.double(rep(0, numReps)))
                        
                        if(p1 > 0){
                            beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                        }
                        if(p1 == 0){
                            beta1.p 		<- NULL
                        }
                        if(p2 > 0){
                            beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                        }
                        if(p2 == 0){
                            beta2.p 		<- NULL
                        }
                        if(p3 > 0){
                            beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                        }
                        if(p3 == 0){
                            beta3.p 		<- NULL
                        }
                        
                        lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                        lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                        lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                        gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                        theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                        mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                        mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                        mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                        sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                        sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                        sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                        K1.p 			<- matrix(mcmcRet$samples_J1, nrow = nStore, byrow = TRUE)
                        K2.p 			<- matrix(mcmcRet$samples_J2, nrow = nStore, byrow = TRUE)
                        K3.p 			<- matrix(mcmcRet$samples_J3, nrow = nStore, byrow = TRUE)
                        s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                        s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                        s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                        
                        if(p1 > 0){
                            accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                        }
                        if(p1 == 0){
                            accept.beta1 	<- NULL
                        }
                        if(p2 > 0){
                            accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                        }
                        if(p2 == 0){
                            accept.beta2 	<- NULL
                        }
                        if(p3 > 0){
                            accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                        }
                        if(p3 == 0){
                            accept.beta3 	<- NULL
                        }
                        accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                        accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                        accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                        accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                        accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                        accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                        accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+7])/sum(as.vector(mcmcRet$moveVec)==11)
                        
                        
                        if(nGam_save > 0 & !is.null(path)){
                            save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                        }
                        
                        if(p1 > 0){
                            covNames1 = colnames(Xmat1)
                        }
                        if(p1 == 0){
                            covNames1 = NULL
                        }
                        
                        if(p2 > 0){
                            covNames2 = colnames(Xmat2)
                        }
                        if(p2 == 0){
                            covNames2 = NULL
                        }
                        if(p3 > 0){
                            covNames3 = colnames(Xmat3)
                        }
                        if(p3 == 0){
                            covNames3 = NULL
                        }
                        
                        
                        ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, theta.p = theta.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3)
                        
                    }	# if(model_h3 == "semi-Markov")
                    
                } # if(hz.type == "PEM")
                
                
                
                chain = chain + 1
                
            }	# while(chain <= nChain)
            
            ret[["setup"]]	<- list(nCov = nCov, hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, nGam_save = nGam_save, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, model = model_h3, nChain = nChain)
            
            
            if(hz.type == "Weibull")
            {
                class(ret) <- c("Bayes", "ID", "Ind", "WB")
            }
            if(hz.type == "PEM")
            {
                class(ret) <- c("Bayes", "ID", "Ind", "PEM")
            }
            
            
            
            return(ret)
            
            #}
            #else{
            #	print("The 'startValues' should be the list of length equal to 'nChain'.")
            #}
        }
        
        
        
        
        ## for cluster-correlated semi-competing risks data
        
        if(!is.null(cluster))
        {
            nChain <- length(startValues)
            
            model_h3 	<- model[1]
            hz.type 	<- model[2]
            re.type 	<- model[3]
            
            Xmat1 <- model.frame(lin.pred[[1]], data=data)
            Xmat2 <- model.frame(lin.pred[[2]], data=data)
            Xmat3 <- model.frame(lin.pred[[3]], data=data)
            
            p1 <- ncol(Xmat1)
            p2 <- ncol(Xmat2)
            p3 <- ncol(Xmat3)
            
            nCov <- c(p1, p2, p3)
            
            
            ##
            survData <- cbind(Y, cluster, Xmat1, Xmat2, Xmat3)
            
            n	<- dim(survData)[1]
            
            J	<- length(unique(survData[,5]))
            
            nj	<- rep(NA, J)
            
            for(i in 1:J){
                nj[i]	<- length(which(survData[,5] == i))
            }
            
            
            #if(class(startValues) == "list" & length(startValues) == nChain){
            
            if(!is.null(path)){
                dir.create(paste(path), recursive = TRUE, showWarnings = FALSE)
            }
            
            
            
            ### setting hyperparameters
            
            if(hz.type == "Weibull" & re.type == "MVN")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab1, hyperP$WB$WB.ab2, hyperP$WB$WB.ab3, hyperP$WB$WB.cd1, hyperP$WB$WB.cd2, hyperP$WB$WB.cd3, hyperP$theta, hyperP$MVN$rho_v, hyperP$MVN$Psi_v))
            }
            if(hz.type == "Weibull" & re.type == "DPM")
            {
                hyperParams <- as.vector(c(hyperP$WB$WB.ab1, hyperP$WB$WB.ab2, hyperP$WB$WB.ab3, hyperP$WB$WB.cd1, hyperP$WB$WB.cd2, hyperP$WB$WB.cd3, hyperP$theta, rep(0, 3), 1, hyperP$DPM$Psi0, hyperP$DPM$rho0, hyperP$DPM$aTau, hyperP$DPM$bTau))
            }
            if(hz.type == "PEM" & re.type == "MVN")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab1, hyperP$PEM$PEM.ab2, hyperP$PEM$PEM.ab3, hyperP$PEM$PEM.alpha1, hyperP$PEM$PEM.alpha2, hyperP$PEM$PEM.alpha3, 1, 1, 1, hyperP$theta, hyperP$MVN$rho_v, hyperP$MVN$Psi_v))
            }
            if(hz.type == "PEM" & re.type == "DPM")
            {
                hyperParams <- as.vector(c(hyperP$PEM$PEM.ab1, hyperP$PEM$PEM.ab2, hyperP$PEM$PEM.ab3, hyperP$PEM$PEM.alpha1, hyperP$PEM$PEM.alpha2, hyperP$PEM$PEM.alpha3, 1, 1, 1, hyperP$theta, rep(0, 3), 1, hyperP$DPM$Psi0, hyperP$DPM$rho0, hyperP$DPM$aTau, hyperP$DPM$bTau))
            }
            
            ### mcmc setting
            
            if(hz.type == "Weibull")
            {
                mcmc <- as.vector(c(mhProp_alpha1_var=mcmcList$tuning$mhProp_alphag_var[1], mhProp_alpha2_var=mcmcList$tuning$mhProp_alphag_var[2], mhProp_alpha3_var=mcmcList$tuning$mhProp_alphag_var[3], mhProp_theta_var=mcmcList$tuning$mhProp_theta_var, mhProp_V1_var=mcmcList$tuning$mhProp_Vg_var[1], mhProp_V2_var=mcmcList$tuning$mhProp_Vg_var[2], mhProp_V3_var=mcmcList$tuning$mhProp_Vg_var[3], nGam_save=mcmcList$storage$nGam_save, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc, storeV=mcmcList$storage$storeV))
            }
            if(hz.type == "PEM")
            {
                mcmcList$tuning$sg_max <- c(1,1,1)
                mcmc <- as.vector(c(C1=mcmcList$tuning$Cg[1], C2=mcmcList$tuning$Cg[2], C3=mcmcList$tuning$Cg[3], delPert1=mcmcList$tuning$delPertg[1], delPert2=mcmcList$tuning$delPertg[2], delPert3=mcmcList$tuning$delPertg[3], rj.scheme = mcmcList$tuning$rj.scheme, K1_max=mcmcList$tuning$Kg_max[1], K2_max=mcmcList$tuning$Kg_max[2], K3_max=mcmcList$tuning$Kg_max[3], s1_max=mcmcList$tuning$sg_max[1], s2_max=mcmcList$tuning$sg_max[2], s3_max=mcmcList$tuning$sg_max[3], mhProp_theta_var=mcmcList$tuning$mhProp_theta_var, mhProp_V1_var=mcmcList$tuning$mhProp_Vg_var[1], mhProp_V2_var=mcmcList$tuning$mhProp_Vg_var[2], mhProp_V3_var=mcmcList$tuning$mhProp_Vg_var[3], nGam_save=mcmcList$storage$nGam_save, numReps=mcmcList$run$numReps, thin=mcmcList$run$thin, burninPerc=mcmcList$run$burninPerc, storeV=mcmcList$storage$storeV))
            }
            
            
            
            chain = 1
            ret <- list()
            
            while(chain <= nChain){
                
                cat("chain: ", chain, "\n")
                nam = paste("chain", chain, sep="")
                
                temp <- startValues[[chain]]
                
                ### setting starting values
                
                if(hz.type == "Weibull" & re.type == "MVN")
                {
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, temp$WB$WB.alpha, temp$WB$WB.kappa, theta=temp$common$theta, gamma=temp$common$gamma.ji, V1=temp$common$V.j1, V2=temp$common$V.j2, V3=temp$common$V.j3, Sigma_V=temp$MVN$MVN.SigmaV))
                }
                if(hz.type == "Weibull" & re.type == "DPM")
                {
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, temp$WB$WB.alpha, temp$WB$WB.kappa, temp$common$theta, gamma=temp$common$gamma.ji, V1=temp$common$V.j1, V2=temp$common$V.j2, V3=temp$common$V.j3, class=temp$DPM$DPM.class, tau=temp$DPM$DPM.tau))
                }
                if(hz.type == "PEM" & re.type == "MVN")
                {
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, theta=temp$common$theta, gamma=temp$common$gamma.ji, V1=temp$common$V.j1, V2=temp$common$V.j2, V3=temp$common$V.j3, Sigma_V=temp$MVN$MVN.SigmaV))
                }
                if(hz.type == "PEM" & re.type == "DPM")
                {
                    startV <- as.vector(c(beta1=temp$common$beta1, beta2=temp$common$beta2, beta3=temp$common$beta3, theta=temp$common$theta, gamma=temp$common$gamma.ji, V1=temp$common$V.j1, V2=temp$common$V.j2, V3=temp$common$V.j3, class=temp$DPM$DPM.class, tau=temp$DPM$DPM.tau))
                }
                
                
                # hz.type = "Weibull"
                
                if(hz.type == "Weibull"){
                    
                    nGam_save   <- mcmc[8]
                    numReps     <- mcmc[9]
                    thin        <- mcmc[10]
                    burninPerc  <- mcmc[11]
                    storeV      <- mcmc[12:14]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    # re.type = "MVN"
                    
                    if(re.type == "MVN"){
                        
                        # model = "Markov"
                        
                        #######################################
                        ############ Weibull-MVN-M ############
                        #######################################
                        
                        if(model_h3 == "Markov"){
                            
                            ###
                            
                            startValues1 <- startV[1:(p1+p2+p3+n+7)]
                            startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmc1 <- mcmc[1:4]
                            mcmc2 <- mcmc[5:7]
                            
                            mcmcNew <- c(mcmc1, 1, mcmc2, n)
                            
                            
                            
                            mcmcRet <- .C("BweibMvnCorScrmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            mcmc		= as.double(mcmcNew),
                            startValues 	= as.double(startV),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_alpha1 	= as.double(rep(0, nStore*1)),
                            samples_alpha2 	= as.double(rep(0, nStore*1)),
                            samples_alpha3 	= as.double(rep(0, nStore*1)),
                            samples_kappa1 	= as.double(rep(0, nStore*1)),
                            samples_kappa2 	= as.double(rep(0, nStore*1)),
                            samples_kappa3 	= as.double(rep(0, nStore*1)),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J+J+J))),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            lpml2     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                            alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                            alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                            kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                            kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                            kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            Sigma_V.p		<- array(as.vector(mcmcRet$samples_Sigma_V), c(3, 3, nStore))
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                            accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                            accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                            accept.V1       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.V2       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+J+1):(p1+p2+p3+7+n+n+J+J)])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.V3       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+J+J+1):(p1+p2+p3+7+n+n+J+J+J)])/sum(as.vector(mcmcRet$moveVec)==15)
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            lpml2.p <- matrix(mcmcRet$lpml2, nrow = nStore, byrow = TRUE)
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                         
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                        } ## end: if Weibull-MVN-M
                        
                        
                        # model = "semi-Markov"
                        
                        #######################################
                        ############ 2 Weibull-MVN-SM #########
                        #######################################
                        
                        if(model_h3 == "semi-Markov"){
                            
                            ###
                            
                            startValues1 <- startV[1:(p1+p2+p3+n+7)]
                            startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmc1 <- mcmc[1:4]
                            mcmc2 <- mcmc[5:7]
                            
                            mcmcNew <- c(mcmc1, 1, mcmc2, n)
                            
                            
                            mcmcRet <- .C("BweibMvnCorScrSMmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            mcmc		= as.double(mcmcNew),
                            startValues 	= as.double(startV),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_alpha1 	= as.double(rep(0, nStore*1)),
                            samples_alpha2 	= as.double(rep(0, nStore*1)),
                            samples_alpha3 	= as.double(rep(0, nStore*1)),
                            samples_kappa1 	= as.double(rep(0, nStore*1)),
                            samples_kappa2 	= as.double(rep(0, nStore*1)),
                            samples_kappa3 	= as.double(rep(0, nStore*1)),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J+J+J))),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            lpml2     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                            alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                            alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                            kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                            kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                            kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            Sigma_V.p		<- array(as.vector(mcmcRet$samples_Sigma_V), c(3, 3, nStore))
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                            accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                            accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                            accept.V1       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.V2       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+J+1):(p1+p2+p3+7+n+n+J+J)])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.V3       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+J+J+1):(p1+p2+p3+7+n+n+J+J+J)])/sum(as.vector(mcmcRet$moveVec)==15)
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            lpml2.p <- matrix(mcmcRet$lpml2, nrow = nStore, byrow = TRUE)
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }  ## end: if Weibull-MVN-SM
                        
                    } ## end: if Weibull-MVN
                    
                    
                    # re.type = "DPM"
                    
                    if(re.type == "DPM"){
                        
                        # model = "Markov"
                        
                        #######################################
                        ############ Weibull-DPM-M ############
                        #######################################
                        
                        if(model_h3 == "Markov"){
                            
                            ##
                            startValues1 <- startV[1:(p1+p2+p3+n+7)]
                            startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmc1 <- mcmc[1:4]
                            mcmc2 <- mcmc[5:7]
                            
                            mcmcNew <- c(mcmc1, 1, mcmc2, n)
                            
                            
                            mcmcRet <- .C("BweibDpCorScrmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            mcmc		= as.double(mcmcNew),
                            startValues 	= as.double(startV),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_alpha1 	= as.double(rep(0, nStore*1)),
                            samples_alpha2 	= as.double(rep(0, nStore*1)),
                            samples_alpha3 	= as.double(rep(0, nStore*1)),
                            samples_kappa1 	= as.double(rep(0, nStore*1)),
                            samples_kappa2 	= as.double(rep(0, nStore*1)),
                            samples_kappa3 	= as.double(rep(0, nStore*1)),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_c		= as.double(rep(0, nStore*J)),
                            samples_mu		= as.double(rep(0, nStore*3*J)),
                            samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                            samples_tau     = as.double(rep(0, nStore*1)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J))),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            lpml2     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                            alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                            alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                            kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                            kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                            kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                            nu2.p           <- matrix(mcmcRet$samples_nu2, nrow = nStore, byrow = TRUE)
                            nu3.p           <- matrix(mcmcRet$samples_nu3, nrow = nStore, byrow = TRUE)
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            c.p            <- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                            mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                            Sigma.p         <- array(as.vector(mcmcRet$samples_Sigma), c(3, 3 * J, nStore))
                            tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                            accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                            accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                            accept.nu2      <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+5)])
                            accept.nu3      <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6)])
                            accept.gamma    <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+1):(p1+p2+p3+6+n)])
                            lastChgProp     <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+n+1)])
                            mhGam_chk       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+n+1+1):(p1+p2+p3+6+n+1+n)])
                            accept.V        <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==13)/3
                            
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            lpml2.p <- matrix(mcmcRet$lpml2, nrow = nStore, byrow = TRUE)
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
      
                            
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, V1.p = V1.p, V2.p = V2.p, V3.p = V3.p, class.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }  ## end: if Weibull-DPM-M
                        
                        
                        # model = "semi-Markov"
                        
                        #######################################
                        ############ Weibull-DPM-SM ############
                        #######################################
                        
                        if(model_h3 == "semi-Markov"){
                            
                            ###
                            startValues1 <- startV[1:(p1+p2+p3+n+7)]
                            startValues2 <- startV[(p1+p2+p3+n+8):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmc1 <- mcmc[1:4]
                            mcmc2 <- mcmc[5:7]
                            
                            mcmcNew <- c(mcmc1, 1, mcmc2, n)
                            
                            mcmcRet <- .C("BweibDpCorScrSMmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            mcmc		= as.double(mcmcNew),
                            startValues 	= as.double(startV),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_alpha1 	= as.double(rep(0, nStore*1)),
                            samples_alpha2 	= as.double(rep(0, nStore*1)),
                            samples_alpha3 	= as.double(rep(0, nStore*1)),
                            samples_kappa1 	= as.double(rep(0, nStore*1)),
                            samples_kappa2 	= as.double(rep(0, nStore*1)),
                            samples_kappa3 	= as.double(rep(0, nStore*1)),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_c		= as.double(rep(0, nStore*J)),
                            samples_mu		= as.double(rep(0, nStore*3*J)),
                            samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                            samples_tau     = as.double(rep(0, nStore*1)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+6+n+1+n+J))),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            lpml2     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            alpha1.p 		<- matrix(mcmcRet$samples_alpha1, nrow = nStore, byrow = TRUE)
                            alpha2.p 		<- matrix(mcmcRet$samples_alpha2, nrow = nStore, byrow = TRUE)
                            alpha3.p 		<- matrix(mcmcRet$samples_alpha3, nrow = nStore, byrow = TRUE)
                            kappa1.p 		<- matrix(mcmcRet$samples_kappa1, nrow = nStore, byrow = TRUE)
                            kappa2.p 		<- matrix(mcmcRet$samples_kappa2, nrow = nStore, byrow = TRUE)
                            kappa3.p 		<- matrix(mcmcRet$samples_kappa3, nrow = nStore, byrow = TRUE)
                            nu2.p           <- matrix(mcmcRet$samples_nu2, nrow = nStore, byrow = TRUE)
                            nu3.p           <- matrix(mcmcRet$samples_nu3, nrow = nStore, byrow = TRUE)
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            c.p            <- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                            mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                            Sigma.p         <- array(as.vector(mcmcRet$samples_Sigma), c(3, 3 * J, nStore))
                            tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.alpha1	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+1)])/sum(as.vector(mcmcRet$moveVec)==4)
                            accept.alpha2	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+2)])/sum(as.vector(mcmcRet$moveVec)==5)
                            accept.alpha3	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+3)])/sum(as.vector(mcmcRet$moveVec)==6)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+4)])/sum(as.vector(mcmcRet$moveVec)==10)
                            accept.nu2      <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+5)])
                            accept.nu3      <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6)])
                            accept.gamma    <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+1):(p1+p2+p3+6+n)])
                            lastChgProp     <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+n+1)])
                            mhGam_chk       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+6+n+1+1):(p1+p2+p3+6+n+1+n)])
                            accept.V        <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7+n+n+1):(p1+p2+p3+7+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==13)/3
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            lpml2.p <- matrix(mcmcRet$lpml2, nrow = nStore, byrow = TRUE)
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                            
                            
                            
                            
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, alpha1.p = alpha1.p, alpha2.p = alpha2.p, alpha3.p = alpha3.p, kappa1.p = kappa1.p, kappa2.p = kappa2.p, kappa3.p = kappa3.p, theta.p = theta.p, V1.p = V1.p, V2.p = V2.p, V3.p = V3.p, class.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.alpha1 = accept.alpha1, accept.alpha2 = accept.alpha2, accept.alpha3 = accept.alpha3, accept.theta = accept.theta, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }  ## end: if Weibull-DPM-SM
                        
                    } ## end: if Weibull-DPM
                    
                } ## end: if Weibull
  
                
                
                # hz.type = "PEM"
                
                if(hz.type == "PEM"){
                    
                    
                    C1 <- mcmc[1]
                    C2 <- mcmc[2]
                    C3 <- mcmc[3]
                    delPert1 <- mcmc[4]
                    delPert2 <- mcmc[5]
                    delPert3 <- mcmc[6]
                    rj.scheme <- mcmc[7]
                    K1_max <- mcmc[8]
                    K2_max <- mcmc[9]
                    K3_max <- mcmc[10]
                    s1_max <- max(temp$PEM$PEM.s1)
                    s2_max <- max(temp$PEM$PEM.s2)
                    s3_max <- max(temp$PEM$PEM.s3)
                    mhProp_theta_var <- mcmc[14]
                    mhProp_V1_var   <- mcmc[15]
                    mhProp_V2_var   <- mcmc[16]
                    mhProp_V3_var   <- mcmc[17]
                    nGam_save   <- mcmc[18]
                    numReps     <- mcmc[19]
                    thin        <- mcmc[20]
                    burninPerc  <- mcmc[21]
                    storeV      <- mcmc[22:24]
                    
                    nStore <- round(numReps/thin*(1-burninPerc))
                    
                    
                    ## recommended when the unit of time is day
                    if(rj.scheme == 1){
                        s_propBI1 <- seq(1, s1_max, 1)
                        s_propBI1 <- s_propBI1[s_propBI1 < s1_max]
                        s_propBI2 <- seq(1, s2_max, 1)
                        s_propBI2 <- s_propBI2[s_propBI2 < s2_max]
                        s_propBI3 <- seq(1, s3_max, 1)
                        s_propBI3 <- s_propBI3[s_propBI3 < s3_max]
                    }
                    ## uniquely ordered failure time
                    if(rj.scheme == 2){
                        s_propBI1 <- sort(unique(survData[survData[,2]==1,1]))
                        s_propBI1 <- s_propBI1[s_propBI1 < s1_max]
                        s_propBI2 <- sort(unique(survData[(survData[,2]==0) & (survData[,4]==1),3]))
                        s_propBI2 <- s_propBI2[s_propBI2 < s2_max]
                        s_propBI3 <- sort(unique(survData[(survData[,2]==1) & (survData[,4]==1),3]))
                        s_propBI3 <- s_propBI3[s_propBI3 < s3_max]
                    }
                    
                    num_s_propBI1=length(s_propBI1)
                    num_s_propBI2=length(s_propBI2)
                    num_s_propBI3=length(s_propBI3)
                    
                    time_lambda1 <- mcmcList$tuning$time_lambda1
                    time_lambda2 <- mcmcList$tuning$time_lambda2
                    time_lambda3 <- mcmcList$tuning$time_lambda3
                    
                    nTime_lambda1 <- length(time_lambda1)
                    nTime_lambda2 <- length(time_lambda2)
                    nTime_lambda3 <- length(time_lambda3)
                    
                    mcmcParams <- c(C1, C2, C3, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2, num_s_propBI3, K1_max, K2_max, K3_max, s1_max, s2_max, s3_max, nTime_lambda1, nTime_lambda2, nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1, time_lambda2, time_lambda3, mhProp_theta_var, mhProp_V1_var, mhProp_V2_var, mhProp_V3_var)
                    
                    s1  		<- temp$PEM$PEM.s1
                    s2			<- temp$PEM$PEM.s2
                    s3			<- temp$PEM$PEM.s3
                    
                    lambda1 <- temp$PEM$PEM.lambda1
                    lambda2 <- temp$PEM$PEM.lambda1
                    lambda3 <- temp$PEM$PEM.lambda3
                    
                    K1=temp$PEM$K1
                    K2=temp$PEM$K2
                    K3=temp$PEM$K3
                    mu_lam1=temp$PEM$PEM.mu_lam[1]
                    mu_lam2=temp$PEM$PEM.mu_lam[2]
                    mu_lam3=temp$PEM$PEM.mu_lam[3]
                    sigSq_lam1=temp$PEM$PEM.sigSq_lam[1]
                    sigSq_lam2=temp$PEM$PEM.sigSq_lam[2]
                    sigSq_lam3=temp$PEM$PEM.sigSq_lam[3]
                    
                    
                    
                    
                    
                    
                    
                    # re.type = "MVN"
                    
                    if(re.type == "MVN"){
                        
                        # model = "Markov"
                        
                        #######################################
                        ############ PEM-MVN-M ############
                        #######################################
                        
                        if(model_h3 == "Markov"){
                            
                            ###
                            K_1 <- K1
                            K_2 <- K2
                            K_3 <- K3
                            
                            theta <- startV[(p1+p2+p3+1)]
                            gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                            
                            startV <- as.vector(c(startV[1:(p1+p2+p3)],
                            K1, K2, K3,
                            mu_lam1, mu_lam2, mu_lam3,
                            sigSq_lam1, sigSq_lam2, sigSq_lam3,
                            theta, gamma,
                            lambda1, lambda2, lambda3, s1, s2, s3,
                            startV[(p1+p2+p3+n+1+1):length(startV)]))
                            
                            
                            startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                            startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                            mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                            
                            mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                            
                            
                            
                            mcmcRet <- .C("BpeMvnCorScrmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            startValues 	= as.double(startV),
                            mcmcParams		= as.double(mcmcParams),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_mu_lam1     = as.double(rep(0, nStore*1)),
                            samples_mu_lam2     = as.double(rep(0, nStore*1)),
                            samples_mu_lam3     = as.double(rep(0, nStore*1)),
                            samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                            samples_K1          = as.double(rep(0, nStore*1)),
                            samples_K2          = as.double(rep(0, nStore*1)),
                            samples_K3          = as.double(rep(0, nStore*1)),
                            samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                            samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                            samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_gamma_last = as.double(rep(0, n)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J+J+J))),
                            lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                            lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                            lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            lpml2     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                            lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                            lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                            
                            mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                            mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                            mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                            sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                            sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                            sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                            
                            K1.p 			<- matrix(mcmcRet$samples_K1, nrow = nStore, byrow = TRUE)
                            K2.p 			<- matrix(mcmcRet$samples_K2, nrow = nStore, byrow = TRUE)
                            K3.p 			<- matrix(mcmcRet$samples_K3, nrow = nStore, byrow = TRUE)
                            s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                            s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                            s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                            
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            Sigma_V.p		<- array(as.vector(mcmcRet$samples_Sigma_V), c(3, 3, nStore))
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                            accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                            accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                            accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7)])/sum(as.vector(mcmcRet$moveVec)==11)
                            accept.V1       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==18)
                            accept.V2       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+J+1):(p1+p2+p3+10+n+n+J+J)])/sum(as.vector(mcmcRet$moveVec)==19)
                            accept.V3       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+J+J+1):(p1+p2+p3+10+n+n+J+J+J)])/sum(as.vector(mcmcRet$moveVec)==20)
                            

                            
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            lpml2.p <- matrix(mcmcRet$lpml2, nrow = nStore, byrow = TRUE)
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                            
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }   ## end: if PEM-MVN-M
                        
                        # model = "semi-Markov"
                        
                        #######################################
                        ############ PEM-MVN-SM ############
                        #######################################
                        
                        if(model_h3 == "semi-Markov"){
                            
                            ###
                            
                            K_1 <- K1
                            K_2 <- K2
                            K_3 <- K3
                            
                            
                            theta <- startV[(p1+p2+p3+1)]
                            gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                            
                            startV <- as.vector(c(startV[1:(p1+p2+p3)],
                            K1, K2, K3,
                            mu_lam1, mu_lam2, mu_lam3,
                            sigSq_lam1, sigSq_lam2, sigSq_lam3,
                            theta, gamma,
                            lambda1, lambda2, lambda3, s1, s2, s3,
                            startV[(p1+p2+p3+n+1+1):length(startV)]))
                            
                            
                            startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                            startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                            mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                            
                            mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                            
                            
                            
                            mcmcRet <- .C("BpeMvnCorScrSMmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            startValues 	= as.double(startV),
                            mcmcParams		= as.double(mcmcParams),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_mu_lam1     = as.double(rep(0, nStore*1)),
                            samples_mu_lam2     = as.double(rep(0, nStore*1)),
                            samples_mu_lam3     = as.double(rep(0, nStore*1)),
                            samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                            samples_K1          = as.double(rep(0, nStore*1)),
                            samples_K2          = as.double(rep(0, nStore*1)),
                            samples_K3          = as.double(rep(0, nStore*1)),
                            samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                            samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                            samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_Sigma_V	= as.double(rep(0, nStore*3*3)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_gamma_last = as.double(rep(0, n)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J+J+J))),
                            lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                            lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                            lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                            lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                            lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                            
                            mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                            mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                            mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                            sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                            sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                            sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                            
                            K1.p 			<- matrix(mcmcRet$samples_K1, nrow = nStore, byrow = TRUE)
                            K2.p 			<- matrix(mcmcRet$samples_K2, nrow = nStore, byrow = TRUE)
                            K3.p 			<- matrix(mcmcRet$samples_K3, nrow = nStore, byrow = TRUE)
                            s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                            s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                            s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                            
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            Sigma_V.p		<- array(as.vector(mcmcRet$samples_Sigma_V), c(3, 3, nStore))
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                            accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                            accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                            accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7)])/sum(as.vector(mcmcRet$moveVec)==11)
                            accept.V1       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==18)
                            accept.V2       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+J+1):(p1+p2+p3+10+n+n+J+J)])/sum(as.vector(mcmcRet$moveVec)==19)
                            accept.V3       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+J+J+1):(p1+p2+p3+10+n+n+J+J+J)])/sum(as.vector(mcmcRet$moveVec)==20)
                            
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, Sigma_V.p = Sigma_V.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V1 = accept.V1, accept.V2 = accept.V2, accept.V3 = accept.V3, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }   ## end: if PEM-MVN-SM
                        
                    }   ## end: if PEM-MVN
                    
                    
                    # re.type = "DPM"
                    
                    if(re.type == "DPM"){
                        
                        # model = "Markov"
                        
                        #######################################
                        ############ PEM-DPM-M ############
                        #######################################
                        
                        if(model_h3 == "Markov"){
                            
                            ###
                            
                            K_1 <- K1
                            K_2 <- K2
                            K_3 <- K3
                            
                            theta <- startV[(p1+p2+p3+1)]
                            gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                            
                            startV <- as.vector(c(startV[1:(p1+p2+p3)],
                            K1, K2, K3,
                            mu_lam1, mu_lam2, mu_lam3,
                            sigSq_lam1, sigSq_lam2, sigSq_lam3,
                            theta, gamma,
                            lambda1, lambda2, lambda3, s1, s2, s3,
                            startV[(p1+p2+p3+n+1+1):length(startV)]))
                            
                            
                            startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                            startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                            mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                            
                            mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                            
                            
                            mcmcRet <- .C("BpeDpCorScrmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            startValues 	= as.double(startV),
                            mcmcParams		= as.double(mcmcParams),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_mu_lam1     = as.double(rep(0, nStore*1)),
                            samples_mu_lam2     = as.double(rep(0, nStore*1)),
                            samples_mu_lam3     = as.double(rep(0, nStore*1)),
                            samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                            samples_K1          = as.double(rep(0, nStore*1)),
                            samples_K2          = as.double(rep(0, nStore*1)),
                            samples_K3          = as.double(rep(0, nStore*1)),
                            samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                            samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                            samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_c		= as.double(rep(0, nStore*J)),
                            samples_mu		= as.double(rep(0, nStore*3*J)),
                            samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                            samples_tau     = as.double(rep(0, nStore*1)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_gamma_last = as.double(rep(0, n)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J))),
                            lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                            lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                            lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                            lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                            lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                            
                            mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                            mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                            mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                            sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                            sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                            sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                            
                            K1.p 			<- matrix(mcmcRet$samples_K1, nrow = nStore, byrow = TRUE)
                            K2.p 			<- matrix(mcmcRet$samples_K2, nrow = nStore, byrow = TRUE)
                            K3.p 			<- matrix(mcmcRet$samples_K3, nrow = nStore, byrow = TRUE)
                            s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                            s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                            s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                            
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            c.p            <- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                            mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                            Sigma.p         <- array(as.vector(mcmcRet$samples_Sigma), c(3, 3 * J, nStore))
                            tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                            accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                            accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                            accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7)])/sum(as.vector(mcmcRet$moveVec)==11)
                            accept.V       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==18)/3
   
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)
                            
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar
                            
                            
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, class.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }   ## end: if PEM-DPM-M
                        
                        
                        
                        # model = "semi-Markov"
                        
                        #######################################
                        ############ PEM-DPM-SM ############
                        #######################################
                        
                        if(model_h3 == "semi-Markov"){
                            
                            ###
   
                            K_1 <- K1
                            K_2 <- K2
                            K_3 <- K3
                            
                            theta <- startV[(p1+p2+p3+1)]
                            gamma <- startV[(p1+p2+p3+2):(p1+p2+p3+n+1)]
                            
                            startV <- as.vector(c(startV[1:(p1+p2+p3)],
                            K1, K2, K3,
                            mu_lam1, mu_lam2, mu_lam3,
                            sigSq_lam1, sigSq_lam2, sigSq_lam3,
                            theta, gamma,
                            lambda1, lambda2, lambda3, s1, s2, s3,
                            startV[(p1+p2+p3+n+1+1):length(startV)]))
                            
                            
                            startValues1 <- startV[1:(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1))]
                            startValues2 <- startV[(p1+p2+p3+10+n+2*(K_1+1)+2*(K_2+1)+2*(K_3+1)+1):length(startV)]
                            startV <- c(startValues1, 1, 1, startValues2)
                            
                            mcmcParams1 <- mcmcParams[1:(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+1)]
                            mcmcParams2 <- mcmcParams[(18+num_s_propBI1+num_s_propBI2+num_s_propBI3+nTime_lambda1+nTime_lambda2+nTime_lambda3+2):(length(mcmcParams)-1)]
                            
                            mcmcParams <- c(mcmcParams1, 1, mcmcParams2, n)
                            
                            mcmcRet <- .C("BpeDpCorScrSMmcmc",
                            survData 		= as.double(as.matrix(survData)),
                            n				= as.integer(n),
                            p1				= as.integer(p1),
                            p2				= as.integer(p2),
                            p3				= as.integer(p3),
                            J				= as.integer(J),
                            nj				= as.double(nj),
                            hyperParams 	= as.double(hyperParams),
                            startValues 	= as.double(startV),
                            mcmcParams		= as.double(mcmcParams),
                            numReps			= as.integer(numReps),
                            thin			= as.integer(thin),
                            burninPerc      = as.double(burninPerc),
                            nGam_save		= as.integer(nGam_save),
                            samples_beta1 	= as.double(rep(0, nStore*p1)),
                            samples_beta2 	= as.double(rep(0, nStore*p2)),
                            samples_beta3 	= as.double(rep(0, nStore*p3)),
                            samples_mu_lam1     = as.double(rep(0, nStore*1)),
                            samples_mu_lam2     = as.double(rep(0, nStore*1)),
                            samples_mu_lam3     = as.double(rep(0, nStore*1)),
                            samples_sigSq_lam1	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam2	= as.double(rep(0, nStore*1)),
                            samples_sigSq_lam3	= as.double(rep(0, nStore*1)),
                            samples_K1          = as.double(rep(0, nStore*1)),
                            samples_K2          = as.double(rep(0, nStore*1)),
                            samples_K3          = as.double(rep(0, nStore*1)),
                            samples_s1          = as.double(rep(0, nStore*(K1_max + 1))),
                            samples_s2          = as.double(rep(0, nStore*(K2_max + 1))),
                            samples_s3          = as.double(rep(0, nStore*(K3_max + 1))),
                            samples_nu2 	= as.double(rep(0, nStore*1)),
                            samples_nu3 	= as.double(rep(0, nStore*1)),
                            samples_theta 	= as.double(rep(0, nStore*1)),
                            samples_V1		= as.double(rep(0, nStore*J)),
                            samples_V2		= as.double(rep(0, nStore*J)),
                            samples_V3		= as.double(rep(0, nStore*J)),
                            samples_c		= as.double(rep(0, nStore*J)),
                            samples_mu		= as.double(rep(0, nStore*3*J)),
                            samples_Sigma	= as.double(rep(0, nStore*3*3*J)),
                            samples_tau     = as.double(rep(0, nStore*1)),
                            samples_gamma 	= as.double(rep(0, nStore*nGam_save)),
                            samples_gamma_last = as.double(rep(0, n)),
                            samples_misc	= as.double(rep(0, (p1+p2+p3+9+n+1+n+J))),
                            lambda1_fin			= as.double(rep(0, nStore*nTime_lambda1)),
                            lambda2_fin			= as.double(rep(0, nStore*nTime_lambda2)),
                            lambda3_fin			= as.double(rep(0, nStore*nTime_lambda3)),
                            gammaP			= as.double(rep(0, n)),
                            dev     			= as.double(rep(0, nStore*1)),
                            invLH = as.double(rep(0, n)),
                            logLH_fin = as.double(0),
                            lpml     			= as.double(rep(0, nStore*1)),
                            moveVec             = as.double(rep(0, numReps)))
                            
                            
                            if(p1 > 0){
                                beta1.p 		<- matrix(mcmcRet$samples_beta1, nrow = nStore, byrow = TRUE)
                            }
                            if(p1 == 0){
                                beta1.p 		<- NULL
                            }
                            if(p2 > 0){
                                beta2.p 		<- matrix(mcmcRet$samples_beta2, nrow = nStore, byrow = TRUE)
                            }
                            if(p2 == 0){
                                beta2.p 		<- NULL
                            }
                            if(p3 > 0){
                                beta3.p 		<- matrix(mcmcRet$samples_beta3, nrow = nStore, byrow = TRUE)
                            }
                            if(p3 == 0){
                                beta3.p 		<- NULL
                            }
                            
                            lambda1.fin 	<- matrix(mcmcRet$lambda1_fin, nrow = nStore, byrow = TRUE)
                            lambda2.fin 	<- matrix(mcmcRet$lambda2_fin, nrow = nStore, byrow = TRUE)
                            lambda3.fin 	<- matrix(mcmcRet$lambda3_fin, nrow = nStore, byrow = TRUE)
                            
                            mu_lam1.p 		<- matrix(mcmcRet$samples_mu_lam1, nrow = nStore, byrow = TRUE)
                            mu_lam2.p 		<- matrix(mcmcRet$samples_mu_lam2, nrow = nStore, byrow = TRUE)
                            mu_lam3.p 		<- matrix(mcmcRet$samples_mu_lam3, nrow = nStore, byrow = TRUE)
                            sigSq_lam1.p 	<- matrix(mcmcRet$samples_sigSq_lam1, nrow = nStore, byrow = TRUE)
                            sigSq_lam2.p 	<- matrix(mcmcRet$samples_sigSq_lam2, nrow = nStore, byrow = TRUE)
                            sigSq_lam3.p 	<- matrix(mcmcRet$samples_sigSq_lam3, nrow = nStore, byrow = TRUE)
                            
                            K1.p 			<- matrix(mcmcRet$samples_K1, nrow = nStore, byrow = TRUE)
                            K2.p 			<- matrix(mcmcRet$samples_K2, nrow = nStore, byrow = TRUE)
                            K3.p 			<- matrix(mcmcRet$samples_K3, nrow = nStore, byrow = TRUE)
                            s1.p 			<- matrix(mcmcRet$samples_s1, nrow = nStore, byrow = TRUE)
                            s2.p 			<- matrix(mcmcRet$samples_s2, nrow = nStore, byrow = TRUE)
                            s3.p 			<- matrix(mcmcRet$samples_s3, nrow = nStore, byrow = TRUE)
                            
                            theta.p 		<- matrix(mcmcRet$samples_theta, nrow = nStore, byrow = TRUE)
                            gamma.p 		<- matrix(mcmcRet$samples_gamma, nrow = nStore, byrow = TRUE)
                            
                            V1.p            <- matrix(mcmcRet$samples_V1, nrow = nStore, byrow = TRUE)
                            V2.p            <- matrix(mcmcRet$samples_V2, nrow = nStore, byrow = TRUE)
                            V3.p            <- matrix(mcmcRet$samples_V3, nrow = nStore, byrow = TRUE)
                            c.p            <- matrix(mcmcRet$samples_c, nrow = nStore, byrow = TRUE)
                            mu.p            <- matrix(mcmcRet$samples_mu, nrow = nStore, byrow = TRUE)
                            Sigma.p         <- array(as.vector(mcmcRet$samples_Sigma), c(3, 3 * J, nStore))
                            tau.p           <- matrix(mcmcRet$samples_tau, nrow = nStore, byrow = TRUE)
                            
                            if(p1 > 0){
                                accept.beta1 	<- as.vector(mcmcRet$samples_misc[1:(p1)])/sum(as.vector(mcmcRet$moveVec)==1)*p1
                            }
                            if(p1 == 0){
                                accept.beta1 	<- NULL
                            }
                            if(p2 > 0){
                                accept.beta2 	<- as.vector(mcmcRet$samples_misc[(p1+1):(p1+p2)])/sum(as.vector(mcmcRet$moveVec)==2)*p2
                            }
                            if(p2 == 0){
                                accept.beta2 	<- NULL
                            }
                            if(p3 > 0){
                                accept.beta3 	<- as.vector(mcmcRet$samples_misc[(p1+p2+1):(p1+p2+p3)])/sum(as.vector(mcmcRet$moveVec)==3)*p3
                            }
                            if(p3 == 0){
                                accept.beta3 	<- NULL
                            }
                            
                            accept.BI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+1])/sum(as.vector(mcmcRet$moveVec)==12)
                            accept.DI1		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+2])/sum(as.vector(mcmcRet$moveVec)==15)
                            accept.BI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+3])/sum(as.vector(mcmcRet$moveVec)==13)
                            accept.DI2		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+4])/sum(as.vector(mcmcRet$moveVec)==16)
                            accept.BI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+5])/sum(as.vector(mcmcRet$moveVec)==14)
                            accept.DI3		<- as.vector(mcmcRet$samples_misc[(p1+p2+p3)+6])/sum(as.vector(mcmcRet$moveVec)==17)
                            accept.theta	<- as.vector(mcmcRet$samples_misc[(p1+p2+p3+7)])/sum(as.vector(mcmcRet$moveVec)==11)
                            accept.V       <- as.vector(mcmcRet$samples_misc[(p1+p2+p3+10+n+n+1):(p1+p2+p3+10+n+n+J)])/sum(as.vector(mcmcRet$moveVec)==18)/3
                            
                            
                            V1summary <- as.matrix(apply(V1.p, 2, summary))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.975))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, quantile, prob = 0.025))
                            V1summary <- rbind(V1summary, apply(V1.p, 2, sd))
                            rownames(V1summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V2summary <- as.matrix(apply(V2.p, 2, summary))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.975))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, quantile, prob = 0.025))
                            V2summary <- rbind(V2summary, apply(V2.p, 2, sd))
                            rownames(V2summary)[7:9] <- c("0.975", "0.025", "sd")
                            
                            V3summary <- as.matrix(apply(V3.p, 2, summary))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.975))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, quantile, prob = 0.025))
                            V3summary <- rbind(V3summary, apply(V3.p, 2, sd))
                            rownames(V3summary)[7:9] <- c("0.975", "0.025", "sd")			
                            
                            if(storeV[1] == TRUE & !is.null(path))
                            {
                                save(V1.p, file = paste(path, "V1Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[2] == TRUE & !is.null(path))
                            {
                                save(V2.p, file = paste(path, "V2Pch", chain, ".RData", sep = ""))
                            }
                            if(storeV[3] == TRUE & !is.null(path))
                            {
                                save(V3.p, file = paste(path, "V3Pch", chain, ".RData", sep = ""))
                            }
                            
                            if(nGam_save > 0 & !is.null(path)){
                                save(gamma.p, file = paste(path, "/gammaPch", chain, ".Rdata", sep = ""))
                            }
                            
                            if(p1 > 0){
                                covNames1 = colnames(Xmat1)
                            }
                            if(p1 == 0){
                                covNames1 = NULL
                            }
                            
                            if(p2 > 0){
                                covNames2 = colnames(Xmat2)
                            }
                            if(p2 == 0){
                                covNames2 = NULL
                            }
                            if(p3 > 0){
                                covNames3 = colnames(Xmat3)
                            }
                            if(p3 == 0){
                                covNames3 = NULL
                            }
                            
                            
                            
                            
                            ### posterior predictive checks ###
                            
                            ## 1. log pseudo-marginal likelihood
                            
                            invLH.p <- matrix(mcmcRet$invLH, nrow = n, byrow = TRUE)		
                            
                            cpo     <- 1/invLH.p
                            
                            LPML <- sum(log(cpo))
                            
                            # or
                            
                            lpml.p <- matrix(mcmcRet$lpml, nrow = nStore, byrow = TRUE)    
                            
                            
                            
                            ## 2. deviance information criterion
                            
                            gamma_mean <- matrix(mcmcRet$gammaP, nrow = n, byrow = TRUE)
                            dev.p 	   <- matrix(mcmcRet$dev, nrow = nStore, byrow = TRUE)	
                            Dbar        <- mean(dev.p)
                            pD          <- Dbar - (-2*mcmcRet$logLH_fin)
                            
                            DIC <- pD + Dbar              
                            
                            
                            
                            
                            ret[[nam]] <- list(beta1.p = beta1.p, beta2.p = beta2.p, beta3.p = beta3.p, lambda1.fin = lambda1.fin, lambda2.fin = lambda2.fin, lambda3.fin = lambda3.fin, mu_lam1.p = mu_lam1.p, mu_lam2.p = mu_lam2.p, mu_lam3.p = mu_lam3.p, sigSq_lam1.p = sigSq_lam1.p, sigSq_lam2.p = sigSq_lam2.p, sigSq_lam3.p = sigSq_lam3.p, K1.p = K1.p, K2.p = K2.p, K3.p = K3.p, s1.p = s1.p, s2.p = s2.p, s3.p = s3.p, theta.p = theta.p, class.p = c.p, mu.p = mu.p, Sigma.p = Sigma.p, tau.p = tau.p, accept.beta1 = accept.beta1, accept.beta2 = accept.beta2, accept.beta3 = accept.beta3, accept.BI1 = accept.BI1, accept.BI2 = accept.BI2, accept.BI3 = accept.BI3, accept.DI1 = accept.DI1, accept.DI2 = accept.DI2, accept.DI3 = accept.DI3, accept.theta = accept.theta, time_lambda1 = time_lambda1, time_lambda2 = time_lambda2, time_lambda3 = time_lambda3, accept.V = accept.V, covNames1 = covNames1, covNames2 = covNames2, covNames3 = covNames3, V1sum = V1summary, V2sum = V2summary, V3sum = V3summary, gamma_mean = gamma_mean)
                            
                            
                        }   ## end: if PEM-DPM-SM
                        
                        
                    }   ## end: if PEM-DPM
                    
                    
                }   ## end: if PEM
                
                
                chain = chain + 1	
                
                
            }	## end: while(chain <= nChain)
            
            
            ret[["setup"]]	<- list(nCov = nCov, hyperParams = hyperParams, startValues = startValues, mcmc = mcmc, nGam_save = nGam_save, numReps = numReps, thin = thin, path = path, burninPerc = burninPerc, hz.type = hz.type, re.type = re.type, model = model_h3, nChain = nChain)
            
            if(hz.type == "Weibull")
            {
                if(re.type == "MVN")
                {
                    class(ret) <- c("Bayes", "ID", "Cor", "WB", "MVN")
                }
                if(re.type == "DPM")
                {
                    class(ret) <- c("Bayes", "ID", "Cor", "WB", "DPM")
                }
            }
            if(hz.type == "PEM")
            {
                if(re.type == "MVN")
                {
                    class(ret) <- c("Bayes", "ID", "Cor", "PEM", "MVN")
                }
                if(re.type == "DPM")
                {
                    class(ret) <- c("Bayes", "ID", "Cor", "PEM", "DPM")
                }
            }


            return(ret)
            
            
            
            #}  ## end: if(class(startValues) == "list" & length(startValues) == nChain)	
            
            
            
            #else{
            #	print("The 'startValues' should be the list of length equal to 'nChain'.")
            #}
        }

        
    }
    else{
    	print(" (numReps * burninPerc) must be divisible by (thin)")
    }
    
    

}# end of function "BayesID"





























