
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  gsFit                   Fits ARMA-GARCH or ARMA-APARCH model  
################################################################################

gsFit <-
function(
    formula = ~ garch(1,1), 
    data,  
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
    include.mean = TRUE, 
    algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"),
    control = NULL,
    tolerance = NULL,
    title = NULL,
    description = NULL)
{  
    # Description:
    #     This functions reads the univariate time series and fits
    #     a ARMA-GARCH/APARCH model with conditional GEV and stable
    #     distribution.
    #     TIME SERIES MODEL:
    #         According to Wurtz et al. (2006) 
    #         xt = mu + a(B)xt + b(B)et,
    #         et = zt * sigmat
    #         zt ~ Dv(0, 1)
    #         sigmat^2 = omega + sum(alphai(e(t-i))^2) + sum(betaj*sigmat-j^2)
    #     REMARKS:
    #     et for ARMA(m,n), n > 1 will be initiated as 0, i.e, e[-n+1:0] = 0.
    #     zt will be initiated as 0.
    #     ht will be initiated as 0.1 (see eq. (22) of Wurtz, 2006).   
    #     VARIABLE NOTATION USED INSIDE THIS FUNCTION:
    #         N: sample size
    #         m,n: ARMA order
    #         p,q, pq = max(p,q): GARCH order. They were called u,v, uv in Wurtz et al.(2006)
    #         mu, a, b: ARMA parameters
    #         omega, alpha, gamma, beta: APARCH vector parameters
    #         garchLLH: Log Likelihood for the ARMA(m,n)-GARCH(p,q) model
    #         sigma: The scale parameter in a pure ARMA(m,n) model with innovations 
    #         D(shape,skew,sigma,location = 0).
    #         llh: The negative of the log-likelihood
    #         parm: ARMA-GARCH parameters concatenated as c(mu,a,b,omega,alpha,gamma,beta)
    #         h: conditional variance. It was referenced as sigmat in model equations
    #         x: data set. It is the time series that will be modelled as ARMA-GARCH
    #         xi: GEV shape parameter. xi > -0.5 (llh) and < 0.5 (finiteness of variance, see Zhao et al. (2011))
    #         AR, MA: if AR (or MA) equals TRUE, then the model has AR (or MA) order equals to zero.
    #         Stable distributions can be specified in all parametrizations: S0, S1 and S2.
    #         Even for stable distribution llh, there's a problem for finding the estimators for
    #         alpha near to 2 and beta > 0. The bad performance on the ARMA-GARCH model was
    #         similar to the single llh estimation of i.i.d samples of stable distribution. 
    #         In our simulations, the param = 2 was chosen to be the ones that performs
    #         better in these situations.
    #     IMPORTANT DETAILS FOR USING THIS FUNCTION:
    #         Some care must be taken when using this function to estimate the parameters
    #         of some models. 
    #         For pure GARCH(p,q) with p,q >= 1 we need to make sure that the 'gamma' variable
    #         is a vector of length 'p' and with all entries equal to zero.
    #         For pure ARCH(p) = GARCH(p,0) we set the variable GARCH equal to TRUE to indicate
    #         that the model has order q = 0. Then, we make q = 1 to estimate the parameters 
    #         of a GARCH(p,1) with beta = 0 and gamma = 0. 
         
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    #   data - vector of data
    #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
    #   include.mean - a logical, should the mean value be estimated ? 
    #   algorithm - the selected algorithm
    #   cond.dist - name of the conditional distribution
    #   title - a character string which allows for a project title
    #   description - a character string which allows for a project description
    
    # Return:
    #   result - An object of class GEVSTABLEGARCH     
      
    # FUNCTION:  

    # Error Treatment on input parameters
    DEBUG = FALSE
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || any(!is.finite(data)))
        stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")     
      
    # Call:
    CALL = match.call()  
  
    # Configuring Tolerance
    if( is.null (tolerance) )
        tolerance = list( TOLG = 1e-8, TOLSTABLE = 1e-2, TOLSTATIONARITY = 1e-3 )
    if( tolerance$TOLSTATIONARITY > 0.05)
        stop("TOLSTATIONARITY can not be > 0.05")
    
    if(is.null(title))
        title = "ARMA-GARCH modelling"
    
    # Getting order model from object formula
    formula.input <- formula
    formula <- .getFormula(formula)
    m <- formula$formula.order[1]
    n <- formula$formula.order[2]    
    p <- formula$formula.order[3]
    q <- formula$formula.order[4]
    APARCH <- formula$isAPARCH
    formula.mean <- ""
    formula.var <- ""
    if(m > 0 || n > 0)
        formula.mean <- formula$formula.mean
    else
        formula.mean <- ""
    if(p > 0)
      formula.var <- formula$formula.var
    else
      formula.var <- ""
    
    # Stop if the algorithm is not supported
    if(algorithm == "sqp.restriction")
     {   
          if(formula.var == "")
              stop("sqp.restriction should only be used with GARCH/APARCH models.")
          if( any ( cond.dist == c("stableS0", "stableS2") ) )
              stop("sqp.restriction algorithm can only be used with \n stable distribution in S1 parametrization, i.e., cond.dist = 'stableS1'.")
    }
    
    # Stop if the user specified a pure arma model
    if(formula.var == "")
          stop("Pure ARMA model not allowed")
    
    
    # Configuring model order
    printRes = TRUE
    ARMAonly <- FALSE
    AR <- FALSE 
    MA <- FALSE 
    GARCH <- FALSE
    if( m == 0) AR <- TRUE
    if( n == 0) MA <- TRUE
    if( q == 0) GARCH <- TRUE
    optim.finished <- FALSE
    if (AR == TRUE)
        m <- 1
    if( MA == TRUE) 
        n <- 1
    if( (p == 0) && (q == 0))
        ARMAonly = TRUE
    if (GARCH == TRUE && !ARMAonly)
        q <- 1
    
    # Initial configurations
    data <- data; 
    N <- length(data)
    out <- NULL # output of the gsFit function
    messages <- NULL # important messages collected during estimation
    out$order <- c(m,n,p,q,include.mean,APARCH)
    TMPvector <- c(if(m != 0 || n != 0) c("arma(",m,",",n,")-"),
                   if(APARCH==FALSE)c("garch(",p,",",q,")"),
                   if(APARCH==TRUE)c("aparch(",p,",",q,")"))
    TMPorder <- paste(TMPvector, collapse="")
    TMPvectorintercept <- c("include.mean:",if(include.mean==TRUE)"TRUE", 
                            if(include.mean==FALSE)"FALSE")
    TMPintercept <- paste(TMPvectorintercept, collapse="")
    out$model <- paste(TMPorder,"##",TMPintercept, collapse="")
    out$cond.dist <- cond.dist
    out$data <- data
    optim.finished <- FALSE
    if(cond.dist == "gat")
        lengthShape = 2
    else
        lengthShape = 1
    
    if(DEBUG)    
        print(c("lengthShape",lengthShape))
    # This function checks if the model is an stationary ARMA process.   
    arCheck <- function(ar) {
        p <- max(which(c(1, -ar) != 0)) - 1
        if (!p) 
          return(TRUE)
        all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
    }


    # BEGIN: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################
    garchLLH = function(parm){
        
        # check if some parameters are NAN
        if(sum(is.nan(parm)) != 0) {return(1e99)}
        
        # Getting parameters from parm vector 
        mu <- parm[1];
        a <- parm[(1+1):(2+m-1)]
        b <- parm[(1+m+1):(2+m+n-1)]
        omega <- parm[1+m+n+1]
        alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
        gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
        beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
        delta <- parm[2+m+n+2*p+q+1] 
        skew <- parm[3+m+n+2*p+q+1]
        shape <- parm[(4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape)]
        
        
        # Configure 'shape' for the 'GAt' distribution
        #if(cond.dist == "gat")
        #    shape <- parm[4+m+n+2*p+q+1,4+m+n+2*p+q+1+1]
        
        # Configuring delta and gamma for Garch estimation
        if( !APARCH) 
        { 
            gm = rep(0,p);
            delta = 2; 
            if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) ) 
                delta = 1
        }
        
        # Setting parameters to accept tapper off MA, AR or GARCH coefficients
        if( AR == TRUE) 
            a <- 0
        if( MA == TRUE) 
            b <- 0
        if( GARCH == TRUE) 
            beta <- 0
        if (include.mean == FALSE)
            mu <- 0
        
        # Avoid being out of parameter space
        cond.general <- FALSE
        cond.normal <- FALSE
        cond.student <- FALSE
        cond.gev <- FALSE
        cond.stable <- FALSE
        parset <- c(omega,alpha,if(!GARCH) beta,delta)
        cond.general <- any(parset < tolerance$TOLG)
        
        if (cond.general || cond.student || cond.gev || cond.stable)
        {
            return(1e99)
        }

        # ARMA stationarity condition check
        if(!arCheck(a))
            return(1e99)

        
        # Filters the Time series to obtain the i.i.d. sequence of 
        # 'innovations' to evaluate the Log-likelihood function
        if(AR == TRUE && MA  == TRUE)
        {
            filteredSeries <- .filterAparch(data = data,p = p,q = q, 
              mu = mu, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[,1]
            hh <- filteredSeries[,2]
        }
        if(AR == FALSE || MA == FALSE)
        {
            filteredArma <- .filterArma(data = data, size.data = N, m = m, n = n, mu = mu, a = a, b = b)
            filteredSeries <- .filterAparch(data = filteredArma,p = p,q = q, 
                              mu = 0, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[,1]
            hh <- filteredSeries[,2]          
          
        }


        # get output Residuals
        if (optim.finished)
        {
            out$residuals <<- as.numeric(z)
            out$sigma.t <<- as.numeric(hh)
            out$h.t <<- as.numeric(hh^delta)           
        }
        
        # Return llh function 
        llh.dens <- .armaGarchDist(z = z, hh = hh, shape = shape, 
                                 skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) || is.na(llh)) 
        {
            llh <- 1e99
        }
        llh
    }
    # END: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################


    # Getting start, lower and upper bound for parameters to perform optimization
    start <- .getStart(data = data,m = m,n = n,p = p,q = q,AR = AR,
                              MA = MA, cond.dist = cond.dist,
                              TOLG = tolerance$TOLG, TOLSTABLE = tolerance$TOLSTABLE)
    if(DEBUG)
    {
        print("start")
        print(start)
    }
    
    # Function that evaluate stationarity conditions to guide parameter estimation.
    garch.stationarity <- function(parm)
    {
      
          # Getting parameters from parm vector 
          # mu <- parm[1];
          # a <- parm[(1+1):(2+m-1)]
          # b <- parm[(1+m+1):(2+m+n-1)]
          omega <- parm[1+m+n+1]
          alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
          gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
          beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
          delta <- parm[2+m+n+2*p+q+1] 
          skew <- parm[3+m+n+2*p+q+1]
          shape <- parm[(4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape)]
          
          model = list(omega = omega, alpha = alpha, gm = gm, beta = beta, 
                       delta = delta, skew = skew, shape = shape)
         
          
          # Return
          .stationarityAparch(model = model, 
                              formula = formula,
                              cond.dist = cond.dist)
    }

    # Optimization procedure using selected algorithms
    modelLLH <- garchLLH


    if (algorithm == "nlminb")
    {              
        fit1 <- nlminb(start[1,], objective = modelLLH,
                       lower = start[2,], upper = start[3,], 
                       control = control)
        out$llh <- fit1$objective
        out$par <- fit1$par
        out$hessian <- fit1$hessian
        out$convergence <- fit1$convergence
    }
    if (algorithm == "nlminb+nm")
    {              
        fit1.partial <- nlminb(start[1,], objective = modelLLH,
                       lower = start[2,], upper = start[3,], 
                       control = control)
        fit1 <- optim(par = fit1.partial$par, fn = modelLLH, 
                      method = "Nelder-Mead", hessian = TRUE)
        out$llh <- fit1$value
        out$par <- fit1$par
        out$hessian <- fit1$hessian
	  out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp")
    {
        fit1 <- solnp(pars = start[1,], fun = modelLLH, 
                    LB = start[2,], UB = start[3,], control = control)
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
        out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp.restriction")
    {
        fit1 <- solnp(pars = start[1,], fun = modelLLH, ineqfun = garch.stationarity, ineqLB = 0,
                    ineqUB = 1-tolerance$TOLSTATIONARITY, LB = start[2,], UB = start[3,], 
                    control = control)
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
        sizeHessian = length(out$hessian[1,])-1
        out$hessian = out$hessian[1:sizeHessian,1:sizeHessian]
        out$convergence <- fit1$convergence
    }
    if ( (any(c("sqp","nlminb")  == algorithm)) && !is.numeric(try(sqrt(diag(solve(out$hessian))), silent = TRUE)))
    {
        # Call the Nelder-Mead method just to calculate the hessian matrix
        # This is due to the fact that the hessian matrix returned
        # by the sqp routine is often difficult to invert, which does not 
        # happen with the hessian returned by the optim routine 
        # using the "Nelder-Mead" algorithm. Note that 
        # we call it with only one interation, so that the algorithm does 
        # not change the parameter values and only calculate the hessian. 
        fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                             method = "Nelder-Mead", hessian = TRUE,
                             control = list(maxit = 1))
        
        # Check if the Nelder-Mead method changed our estimated parameters
        if(sum(abs(fit1.partial$par-fit1$par)) == 0)
        {   
            out$hessian = fit1.partial$hessian
            if(DEBUG)
                print(abs(fit1.partial$par-fit1$par))
        }
    }
    if ( (algorithm == "sqp.restriction") )
    {
      # Call the Nelder-Mead method just to calculate the hessian matrix
      # This is due to the fact that the hessian matrix returned
      # by the sqp routine is often difficult to invert, which does not 
      # happen with the hessian returned by the optim routine 
      # using the "Nelder-Mead" algorithm. Note that 
      # we call it with only one interation, so that the algorithm does 
      # not change the parameter values and only calculate the hessian. 
      fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                            method = "Nelder-Mead", hessian = TRUE,
                            control = list(maxit = 1))
      
      # Check if the Nelder-Mead method changed our estimated parameters
      if(sum(abs(fit1.partial$par-fit1$par)) == 0)
      {   
        out$hessian = fit1.partial$hessian
        if(DEBUG)
          print(abs(fit1.partial$par-fit1$par))
      }
    }
    
    if(DEBUG)
      print(fit1)

    # Configuring the convergence variable
    if(out$llh == 1e99)
       out$convergence = 1

    # Organizing the output of the program
    optim.finished = TRUE
    
    # Call garchLLH function to update the values of the ARMA residuals 
    # and the GARCH/APARCH volatility.
    modelLLH(out$par)
    
    # Creating index to create a vector with the estimated parameters.
    if(!ARMAonly)
    {   
        outindex <- c(if(include.mean) 1, 
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","skstd","stableS0","stableS1","stableS2","gat")  == cond.dist)) (3+m+n+2*p+q+1),
                    if(any(c("std","gev","stableS0","stableS1","stableS2","sstd","skstd","ged","gat")  == cond.dist)) 
                      (4+m+n+2*p+q+1):(4+m+n+2*p+q+lengthShape))
    } else {
        outindex <- c(if(include.mean) 1, 
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","skstd","stableS0", "stableS1", "stableS2","gat")  == cond.dist)) (1+m+n+2*p+q+1),
                    if(any(c("std","gev","stableS0", "stableS1", "stableS2","sstd","skstd","ged","gat")  == cond.dist)) 
                      (2+m+n+2*p+q+1):(2+m+n+2*p+q+lengthShape),
                    length(out$par))  
    }
    
    outnames <- c(if(include.mean) "mu", 
                  if( AR == FALSE) paste("ar", 1:m, sep = ""),
                  if(MA == FALSE) paste("ma", 1:n, sep = ""),
                  if(!ARMAonly) "omega",
                  if(!ARMAonly) paste("alpha", 1:p, sep = ""),
                  if(APARCH) paste("gamma", 1:p, sep = ""),
                  if(!GARCH) paste("beta", 1:q, sep = ""),
                  if(APARCH) "delta",
                  if(any(c("sstd","stableS0", "stableS1", "stableS2","gat","skstd")  == cond.dist)) "skew",
                  if(any(c("std","gev","stableS0", "stableS1", "stableS2","sstd","ged","gat","skstd")  == cond.dist)) 
                    paste("shape", 1:lengthShape, sep = ""),
                  if(ARMAonly) "sigma")
    
    if(DEBUG)
    {
        print(c("out",out))
        print(c("outindex",outindex))
    }
      
    out$par <- out$par[outindex]
    names(out$par) <- outnames
    out$hessian <- out$hessian[outindex,outindex]
    nParam <- length(out$par)
    out$aic  = 2*out$llh + 2*nParam # out$llh is the negative of the log-likelihood
    out$aicc = 2*out$llh + 2*(nParam+1)*N/(N - nParam - 2)
    out$bic =  2*out$llh + nParam*log(N)
    out$ics = c(out$aic,out$bic,out$aicc)
    names(out$ics) <- c("AIC","BIC","AICc")
    
    # Print Summary
    if ( printRes ) 
    {
        solveHessianFailed = FALSE
        out$se.coef <- 0
        out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
        if(!is.numeric(out$se.coef))
        {
            solveHessianFailed = TRUE
            messages$solving.hessian.matrix = "Error inverting Hessian Matrix"
            out$matcoef <- cbind(out$par, rep(NA,length(out$par)), rep(NA,length(out$par)), rep(NA,length(out$par)))
            dimnames(out$matcoef) = dimnames(out$matcoef) = 
            list(names(out$par), c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
        }
        else
        {
            out$tval <- try(out$par/out$se.coef, silent = TRUE)
            out$matcoef = cbind(out$par, if(is.numeric(out$se.coef)) out$se.coef, if(is.numeric(out$tval)) out$tval, 
                                if(is.numeric(out$tval)) 2*(1-pnorm(abs(out$tval))))
            dimnames(out$matcoef) = list(names(out$tval), 
                                    c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
        }
        cat("\nFinal Estimate of the Negative LLH:\n")
        cat("-LLH:",out$llh)

        if(out$convergence == 0)
            messages$optimization.algorithm = "Algorithm achieved convergence"
        else 
            messages$optimization.algorithm = "Algorithm did not achieved convergence"

        cat("\nCoefficient(s):\n")
        printCoefmat(round(out$matcoef,digits=6), digits = 6, signif.stars = TRUE)
    }

    out$order <- c(formula$formula.order[1],formula$formula.order[2],
                   formula$formula.order[3],formula$formula.order[4])
    names(out$order) <- c("m","n","p","q")
    fit <- list(par = out$par, llh = out$llh, hessian = out$hessian, ics = out$ics,
                order = out$order, cond.dist = cond.dist, se.coef = out$se.coef,
                tval = out$tval, matcoef = out$matcoef)

    # creating the output object as in fGarch package...
    new("GEVSTABLEGARCH", call = as.call(match.call()), formula = formula.input, 
    method = "Max Log-Likelihood Estimation", 
    convergence = out$convergence,
    messages = messages,
    data = data, fit = fit, residuals = out$residuals,
    h.t = out$h.t, sigma.t = as.vector(out$sigma.t), title = as.character(title), 
    description = as.character(description))
}

# -----------------------------------------------------------------------------




GSgarch.dstable <<- function(x,alpha = 1.5, beta = 0, gamma = 1, 
                             delta = 0, param = 1)
{
  return(stabledist::dstable(x, alpha, beta, gamma, 
                             delta, pm = param))
}
if(getOption('.stableIsLoaded', default = FALSE) == TRUE)
{
  GSgarch.dstable <<- function(x,alpha = 1.5, beta = 0, gamma = 1, 
                               delta = 0, param = 1)
  {
    return(stable::dstable.quick(x, alpha, beta, gamma, 
                                 delta, param))
  }
}



# -----------------------------------------------------------------------------
