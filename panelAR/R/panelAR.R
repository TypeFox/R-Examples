### Main panelAR function
### Author: Konstantin Kashin
### August 25, 2013

panelAR <- function(formula, data, panelVar, timeVar, autoCorr = c("ar1", 
    "none", "psar1"), panelCorrMethod = c("none","phet","pcse","pwls","parks"), rhotype ="breg", bound.rho = FALSE, rho.na.rm = FALSE, panel.weight = c("t-1", "t"), dof.correction = FALSE, complete.case = FALSE, seq.times = FALSE,  singular.ok=TRUE) {
    # save environment
    env.base <- environment()
    
    # save call
    call <- match.call()
    
    # check for presence of a formula and data object
    ind <- match(c("formula","data","timeVar"),names(call),nomatch=0)
    if(ind[1]==0){
    	stop("A formula argument is required.",call.=FALSE)
    }
    if(ind[2]==0){
    	stop("A data argument is required.",call.=FALSE)
    }
    if(ind[3]==0){
    	stop("You must specify time ID using timeVar.",call.=FALSE)
    }
    
    # identify autocorrelation and panel structure
    autoCorr.method <- match.arg(autoCorr)
    panelCorr.method <- match.arg(panelCorrMethod)
    
    # identify 2nd stage method
    pMethod <- switch(panelCorr.method, "none"="OLS", "phet"="OLS", "pcse"="OLS", "pwls"="GLS", "parks"="GLS")
        
    # check rhotype
    rhotype <- match.arg(rhotype,c("breg","freg","dw","theil-nagar","scorr","theil"))
    
    # check panel.weight & rho.na.rm arguments
    panel.weight <- match.arg(panel.weight) 
    if(!is.logical(rho.na.rm)){
    	stop("rho.na.rm should be logical argument.")
    }
    
    # make sure that timevar is in data object; check for NA values
   	if (!(timeVar %in% colnames(data))) {
        stop("Please make sure that timeVar is in the data object.",call.=FALSE)
    }
    # extract time vector and check for NA values
    time.vec <- data[, timeVar]
    if (any(is.na(time.vec))) {
       stop("You cannot have NA values for the time variable.",call.=FALSE)
    }
    if (!is.integer(time.vec) & !is.numeric(time.vec)) {
       stop("The time variable must be defined as an integer.",call.=FALSE)
    }

    # extract panel variable and check if is.null
    if (!is.null(panelVar)){
    	if (!(panelVar %in% colnames(data))) {
        stop("Please make sure that panelVar is in the data object.",call.=FALSE)
    	} else{
    	 panel.vec <- as.character(as.vector(data[, panelVar]))
    	 if (any(is.na(panel.vec))) {
    		stop("You cannot have NA values for the panel ID variable.",call.=FALSE)
    	 }
    	}
    } else{   
    	panel.vec <- rep("1",nrow(data))
    	if(panelCorr.method!="none") warning("Without specifying a panel ID variable, data is assumed to come from a single panel and variance is assumed to be homoskedastic within that panel. Panel heteroskedasticity and/or correlation is ignored.",call.=FALSE)
    }
    
    # check for duplicate times in each panel
    if (!is.null(timeVar)) {
        if (any(by(time.vec, panel.vec, function(x) any(table(x) > 1)))) {
           stop("You must specify unique times for each observation in a given panel.",call.=FALSE)
        }
    }
    
    # sort dataframe
    order.index <- order(panel.vec, time.vec)
    data <- data[order.index, ]
    panel.vec <- panel.vec[order.index]
    time.vec <- time.vec[order.index]
    
    # run OLS and extract results
    lm.out <- lm(formula = formula, data = data,singular.ok=singular.ok)
    
    # extract terms
    mterms <- lm.out$terms
    # aliased coefficients
    aliased <- is.na(coef(lm.out))
    # get model matrix, frame, and response from OLS regression
    X <- model.matrix(lm.out)[,!aliased]
    mf <- model.frame(lm.out)
    y <- model.response(mf)
    # create yX matrix
    yX <- cbind(y, X)
    # extract initial residuals from OLS regression
    original.e <- residuals(lm.out)
    # save variable names
    var.names <- colnames(X)
    # number of observations in regression
    N <- length(y)
    
	# get rank and residual dof (rdf)
	rank <- lm.out$rank
	rdf <- N-rank
	    
    # lm automatically row-deletes missing data 
    # so we need to do same for panelVar and timeVar
    # store na.action
    obs.dropped <- lm.out$na.action
    if (!is.null(obs.dropped)) {
        data <- data[-obs.dropped, ]
        panel.vec <- panel.vec[-obs.dropped]
        time.vec <- time.vec[-obs.dropped]
    }
    
    # generate time sequentially if seq.times
    if (seq.times) {
 		time.vec <- as.vector(unlist(by(data, panel.vec, function(x) 1:nrow(x))))
        data[, timeVar] <- time.vec # not sure if we need this
    }
    
    # sort units and times
    units <- sort(unique(panel.vec))
    times <- sort(unique(time.vec))
    # number of units and times
    N.units <- length(units)
    N.times <- length(times)
    # average number of units per panel
    N.avgperpanel <- N/N.units
    
    # stop if just 1 time period
    if(N.times<2){
    	stop("More than two time periods required.",call.=FALSE)
    }
    
    # check on Parks method
    if(panelCorr.method=="parks" & (N.units > N.times)){
    	stop("Cannot estimate Parks-Kmenta method because of singularity.")
    }
    
    # expected number of observations if panels balanced
    NT <- N.units * N.times
    # check if balanced
    balanced <- ifelse(N == NT, TRUE, FALSE)
    
    ### create observations matrix: rows are units and columns are times
    # dimensions: N_p x T
    # if cell is TRUE, panel i at time t is observed
    obs.mat <- reshape(cbind(data[, c(panelVar, timeVar)], TRUE), timevar = timeVar, 
    idvar = panelVar, direction = "wide", new.row.names = units)[,-1]
    col.order <- order(as.integer(gsub("TRUE.", "", colnames(obs.mat))))
    obs.mat <- obs.mat[, col.order]
    colnames(obs.mat) <- times
    obs.mat[is.na(obs.mat)] <- FALSE
    obs.mat <- as.matrix(obs.mat)
    
    ### runs analysis: generate number of runs by panel
    # output is vector of counts of number of runs per panel
    # length of vector is N_p
    N.runs.panel <- apply(obs.mat, MARGIN = 1, function(x) sum(rle(x)$values == TRUE))
    if (any(N.runs.panel > 1)) {
        message(paste("The following units have non-consecutive observations. Use runs.analysis() on output for additional details: ", paste(names(N.runs.panel[N.runs.panel>1]),collapse=", "),".",sep=""))
    }
    
    ### reshape residuals into T x N_p matrix
    # matrix will have NAs if unbalanced design
    e.mat <- matrix(NA, nrow = ncol(obs.mat), ncol = nrow(obs.mat))
    e.mat[t(obs.mat) == TRUE] <- original.e
   
   	### Run Prais-Winsten correction
    pw.output <- prais.correct(method = autoCorr.method,env.base=env.base)

    ### Extract residuals from Prais-Winsten regression
    transformed.resids <- pw.output$pw.lm$residuals
    ### Extract model matrix and response from Prais-Winsten regression
    model.mat.pw <- model.matrix(pw.output$pw.lm)
        
    ### In cases where rhos not bounded to [-1,1], we lose first observation
    ### so need to remake observation matrix
    # observation matrix is now T x N_p
    obs.mat.pw <- t(obs.mat)
    
    # create new observation matrix and update panel and time ID vectors 
    if (!is.null(pw.output$pw.lm$na.action)) {
        panel.vec.pw <- panel.vec[-pw.output$pw.lm$na.action]
        time.vec.pw <- time.vec[-pw.output$pw.lm$na.action]
        obs.mat.pw[obs.mat.pw == TRUE][pw.output$pw.lm$na.action] <- FALSE
    } else {
        panel.vec.pw <- panel.vec
        time.vec.pw <- time.vec
    }
    
    ### Panel homoskedasticity
    if (panelCorr.method == "none") {
    		# calculate residual variance
        sigma <- mean(transformed.resids^2)
        # N_p x N_p matrix of panel covariances
        Sigma <- diag(sigma, N.units)
        # matrix of residual covariances
        Omega <- diag(sigma, nrow(model.mat.pw))
        # number of panel covariances calculated
        N.cov <- 1
        # run estimation
        res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
      ### Panel heteroskedasticity
      } else if (panelCorr.method == "phet"|panelCorr.method == "pwls") {
      	# N_p x N_p matrix of panel covariances
      	sigma.vec <- as.vector(by(transformed.resids, panel.vec.pw, function(x) mean(x^2)))
        if(length(sigma.vec)>1){
        	Sigma <- diag(sigma.vec)
        } else{
        	Sigma <- sigma.vec
        	}
        # matrix of residual covariances
        Omega <- diag(rep(sigma.vec, times = as.integer(table(panel.vec.pw))))
        # number of panel covariances calculated
        N.cov <- length(unique(panel.vec.pw))
        # run estimation
        res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
    } else {
      
      ### Correlated panels / PCSE
        # if balanced panel design
        if (balanced) {
            # set up new T x N_p matrix for residuals
            E <- matrix(transformed.resids, nrow = N.times, ncol = N.units, byrow = FALSE)
            # E'E
            E.E <- crossprod(E)
            # N_p x N_p matrix which 
            # gives number of overlapping time observations between each panel
            weight.mat <- crossprod(obs.mat.pw)
        } else {
            
            # set up new T x N_p matrix for residuals
            E <- obs.mat.pw
            # set missing residuals to 0
            E[E == TRUE] <- transformed.resids
            E[E == FALSE] <- 0
            
            # complete case calculation of E'E and weight.mat
            if (complete.case) {
            		# average times per panel
            	    # select only those time periods with complete cases
                I.com.case <- apply(obs.mat.pw, MARGIN = 1, function(x) prod(x) == 1)
                if (!any(I.com.case)) {
                  stop("Unable to compute correlated SEs / PCSEs because there are no time periods in common across all units. Instead, consider setting complete.case=FALSE.",call.=FALSE)
                } else {
                  if(sum(I.com.case)<(0.5*N.avgperpanel)){
                  	warning(paste("The number of time periods used for the calculation of correlated SEs / PCSEs (",as.character(sum(I.com.case)),") is less than half the average number of time periods per panel (",as.character(round(N.avgperpanel,digits=2)),"). Consider setting complete.case=FALSE.",sep=""),call.=FALSE)
                  }
                  E[!I.com.case, ] <- 0
                  E.E <- crossprod(E)
                  weight.mat <- matrix(data = sum(I.com.case), nrow = nrow(E.E), ncol = nrow(E.E))
                }
            } else {
            	# pairwise option
                E.E <- crossprod(E)
                weight.mat <- crossprod(obs.mat.pw)
            }
        }
        
        # N_p x N_p matrix of panel covariances
        Sigma <- E.E/weight.mat
        # Number of panel covariances calculated
        N.cov <- length(Sigma[lower.tri(Sigma, diag = TRUE)])
        # assume covariance is 0 between panels that do not overlap in time periods
        Sigma <- replace(Sigma, is.na(Sigma), 0)
	
		# matrix of residual covariances
        Omega <- kronecker(Sigma, diag(1, N.times))
        if (!balanced) {
            Omega <- Omega[as.vector(obs.mat.pw), as.vector(obs.mat.pw)]
        }
        
        # run estimation
        res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
    }
	# end of correlated SE calculations    

	# clean up coef vector and var-covar matrix
    coef <- as.vector(res$coef)
    vcov <- res$vcov
    colnames(vcov) <- rownames(vcov) <- names(coef) <- var.names
    
    # calculate yhat
    yhat <- as.vector(X %*% coef)
    names(yhat) <- row.names(X)
    
    # calculate residuals as y-yhat
    resids <- y - yhat
    
    # dof correction to vcov
    if (dof.correction) {
        vcov <- vcov * (N/(N - rank))
    }
    
    # create panelStructure list
    panelStructure <- list(obs.mat = obs.mat, rho = pw.output$pw.rho, Sigma = Sigma, N.cov=N.cov)
    
    if(autoCorr.method=="psar1"){
    names(panelStructure$rho) <- rownames(obs.mat)
    	}
    
    # R2 (if ols)
    if(pMethod=="OLS"){
     transform.y.vec <- model.response(model.frame(pw.output$pw.lm))
   	 r2 <- 1 - sum(transformed.resids^2)/sum((transform.y.vec-mean(transform.y.vec))^2)
    } else{
    		r2 <- NULL
    	}
    
    # insert panelVar and timeVar into model frame for output
    mf[,panelVar] <- panel.vec
    mf[,timeVar] <- time.vec
    
    # create output list
    fit <- list(coefficients = coef, residuals = resids,  fitted.values = yhat, rank = rank,
    				df.residual = rdf,call = call,terms = mterms, model = mf,
    				aliased=aliased, na.action = obs.dropped, vcov=vcov, 
    				r2=r2, panelStructure = panelStructure)
    class(fit) <- "panelAR"
    return(fit)
}  # end of panelAR fxn