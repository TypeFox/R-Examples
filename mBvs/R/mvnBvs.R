
mvnBvs  <-  function(Y,
                    lin.pred,
					data,
                    model="unstructured",
					hyperParams,
					startValues,
					mcmcParams)
{
    numReps     <- mcmcParams$run$numReps
    thin        <- mcmcParams$run$thin
    burninPerc  <- mcmcParams$run$burninPerc
	nStore      <- numReps/thin * (1 - burninPerc)
    
    if((numReps / thin * burninPerc) %% 1 == 0)
    {
        nChain <- length(startValues)
        
        chain = 1
        ret <- list()
        
        while(chain <= nChain)
        {
            cat("chain: ", chain, "\n")
            nam = paste("chain", chain, sep="")
            temp <- startValues[[chain]]
            
            Xmat.vs <- model.frame(lin.pred[[1]], data=data)
            Xmat.adj <- model.frame(lin.pred[[2]], data=data)
            
            ###
            n		<- dim(Y)[1]
            p		<- dim(Y)[2]
            q_adj	<- ncol(Xmat.adj)
            q		<- ncol(Xmat.vs) + q_adj
            
            Data <- cbind(Y, Xmat.vs, Xmat.adj)
            
            if(q_adj > 0){
                covNames = c(colnames(Xmat.vs), colnames(Xmat.adj))
            }
            if(q_adj == 0){
                covNames = colnames(Xmat.vs)
            }
            
            
            B <- matrix(c(temp[(p+1):(q*p + p)]), ncol = p, byrow = F)
            gamma <- matrix(c(temp[(q*p + p+1):(q*p + p + q*p)]), ncol = p, byrow = F)
            
            if(q_adj > 0)
            {
                for(i in 1:q_adj){
                    gamma[q - i + 1,] <- 1
                }
            }
            
            rwBetaVar <- mcmcParams$tuning$mhProp_beta_var
            
            if(model == "factor-analytic")
            {
                startV <- c(temp[p+p*q+p*q+1], temp[1:p], temp[(p+p*q+p*q+1+1):(p+p*q+p*q+1+p)])
                
                mcmcP <- c(mcmcParams$tuning$mhProp_lambda_var, 1, 1, 5)
                hyperP <- c(hyperParams$eta, hyperParams$beta0[p+1], hyperParams$FA$lambda[p+1], hyperParams$FA$sigmaSq[1],  hyperParams$FA$sigmaSq[2], hyperParams$beta0[1:p], hyperParams$FA$lambda[1:p], hyperParams$v, hyperParams$omega)
                
                mcmc <- .C("MBVSfamcmc",
                Data            	= as.double(as.matrix(Data)),
                n                	= as.integer(n),
                p                	= as.integer(p),
                q                	= as.integer(q),
                q_adj               = as.integer(q_adj),
                hyperParams         = as.double(hyperP),
                startValues         = as.double(startV),
                startB              = as.double(B),
                startGamma          = as.double(gamma),
                mcmcParams          = as.double(mcmcP),
                rwBetaVar			= as.double(rwBetaVar),
                numReps             = as.integer(numReps),
                thin                = as.integer(thin),
                burninPerc          = as.double(burninPerc),
                samples_beta0 		= as.double(rep(0, nStore*p)),
                samples_lambda      = as.double(rep(0, nStore*p)),
                samples_sigSq       = as.double(rep(0, nStore*1)),
                samples_B           = as.double(rep(0, nStore*p*q)),
                samples_gamma       = as.double(rep(0, nStore*p*q)),
                samples_misc        = as.double(rep(0, p*q+p+1)))
                
                beta0.p 		<- matrix(mcmc$samples_beta0, nrow = nStore, byrow = TRUE)
                lambda.p 		<- matrix(mcmc$samples_lambda, nrow = nStore, byrow = TRUE)
                sigSq.p 		<- matrix(mcmc$samples_sigSq, nrow = nStore, byrow = TRUE)
                B.p				<- array(as.vector(mcmc$samples_B), c(q, p, nStore))
                gamma.p			<- array(as.vector(mcmc$samples_gamma), c(q, p, nStore))
                
                accept.B		<- matrix(as.vector(mcmc$samples_misc[1:(p*q)]), nrow = q, byrow = FALSE)
                accept.lambda 	<- as.vector(mcmc$samples_misc[(p*q+1):(p*q+p)])
                
                ret[[nam]] <- list(B.p = B.p, gamma.p = gamma.p, beta0.p = beta0.p, lambda.p = lambda.p, sigSq.p = sigSq.p, accept.B = accept.B, accept.lambda = accept.lambda, covNames = covNames)
                
            }else if(model == "unstructured")
            {
                startV <- temp
                Sigma <- matrix(c(temp[(p+q*p+q*p+1):(p+q*p+q*p+p*p)]), ncol = p, byrow = F)
                
                mcmcP <- c(mcmcParams$tuning$mhPsi_prop, mcmcParams$tuning$mhrho_prop, mcmcParams$tuning$mhProp_beta_var)
                hyperP <- c(hyperParams$eta, hyperParams$beta0[p+1], hyperParams$beta0[1:p], hyperParams$v, hyperParams$omega, hyperParams$US$US.Sigma[1], hyperParams$US$US.Sigma[2:(p*p+1)])
                
                mcmc <- .C("MBVSusmcmc",
                Data            	= as.double(as.matrix(Data)),
                n                	= as.integer(n),
                p                	= as.integer(p),
                q                	= as.integer(q),
                q_adj               = as.integer(q_adj),
                hyperParams         = as.double(hyperP),
                startValues         = as.double(startV),
                startB              = as.double(B),
                startGamma          = as.double(gamma),
                startSigma          = as.double(Sigma),
                mcmcParams          = as.double(mcmcP),
                rwBetaVar			= as.double(rwBetaVar),
                numReps             = as.integer(numReps),
                thin                = as.integer(thin),
                burninPerc          = as.double(burninPerc),
                samples_beta0 		= as.double(rep(0, nStore*p)),
                samples_B           = as.double(rep(0, nStore*p*q)),
                samples_gamma       = as.double(rep(0, nStore*p*q)),
                samples_Sigma       = as.double(rep(0, nStore*p*p)),
                samples_misc        = as.double(rep(0, p*q+1)))
                
                beta0.p 		<- matrix(mcmc$samples_beta0, nrow = nStore, byrow = TRUE)
                B.p				<- array(as.vector(mcmc$samples_B), c(q, p, nStore))
                gamma.p			<- array(as.vector(mcmc$samples_gamma), c(q, p, nStore))
                Sigma.p			<- array(as.vector(mcmc$samples_Sigma), c(p, p, nStore))
                
                accept.B		<- matrix(as.vector(mcmc$samples_misc[1:(p*q)]), nrow = q, byrow = FALSE)
                accept.Sigma	<- as.vector(mcmc$samples_misc[p*q+1])
                
                ret[[nam]] <- list(B.p = B.p, gamma.p = gamma.p, Sigma.p = Sigma.p, beta0.p = beta0.p, accept.B = accept.B, accept.Sigma = accept.Sigma, covNames = covNames)
            }

            chain = chain + 1
        }
        
        ret[["setup"]]	<- list(hyperParams = hyperParams, startValues = startValues, mcmcParams = mcmcParams, numReps = numReps, thin = thin, burninPerc = burninPerc, model = model, nChain = nChain)
        
        if(model == "factor-analytic")
        {
            class(ret) <- c("mvnBvs", "factor-analytic")
        }
        if(model == "unstructured")
        {
            class(ret) <- c("mvnBvs", "unstructured")
        }
        return(ret)
    }
    else
    {
        print(" (numReps * burninPerc) must be divisible by (thin)")
    }



}
