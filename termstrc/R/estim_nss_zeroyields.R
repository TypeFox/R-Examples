#########################################################################
### Nelson/Siegel-type yield curve estimation method for 'zeroyields' ###
#########################################################################

estim_nss.zeroyields <- function (dataset,
                                  method = "ns",
                                  lambda = 0.0609*12,
                                  tauconstr = NULL,
                                  optimtype = "firstglobal",
                                  constrOptimOptions = list(control = list(), outer.iterations = 200, outer.eps = 1e-04),...)
  {
    obj <- dataset
    optparam <- matrix(nrow = nrow(obj$yields), ncol = length(get_paramnames(method))) 
    opt_result <- list()
    spsearch <- list()
    
      ## default tau constraints (if not specified by user)
      if (is.null(tauconstr) && method != "dl"){
        tauconstr <- c(min(obj$maturities), max(obj$maturities), 0.1, 0.5)
        if (method == "asv") {tauconstr[4] = 0}
        if (method == "ns") {tauconstr = tauconstr[1:3]}
        print("The following constraints are used for the tau parameters:")
        print(tauconstr)
      }
      
      objfct <- get_objfct(method)
      grad_objfct <- get_grad_objfct(method)
      
      sp_search <- findstartparamyields(obj$yields[1,],obj$maturities, method, tauconstr)
      startparam <- sp_search$startparam
      constraints <- get_constraints(method, tauconstr)

      ## Estimation loop
      
      for (i in 1:nrow(obj$yields)){  
        yields <- obj$yields[i,]
        
        if (i==1) {
          beta <- startparam
          spsearch[[i]] <- sp_search
        }
        
        if(i>1 && optimtype == "allglobal"){
          sp_search <- findstartparamyields(yields,obj$maturities, method, tauconstr)
          beta <- sp_search$startparam
          spsearch[[i]] <- sp_search
        }
        if(i>1 && optimtype == "firstglobal"){
          beta <- optparam[i-1,]
        }
        
        opt_result[[i]] <- estimateyieldcurve(method, yields, obj$maturities, beta, lambda, objfct,
                                              grad_objfct, constraints, constrOptimOptions)
        optparam[i,] <- opt_result[[i]]$par
      }
    
    colnames(optparam) <- get_paramnames(method)
    
    
    yhat <- t(apply(optparam,1, function(x) spotrates(method,x,obj$maturities,lambda)))
    
    result <- list(optparam = optparam, opt_result = opt_result, method = method,
                   maturities = obj$maturities, dates = obj$dates, spsearch = spsearch,
                   yields = obj$yields,yhat = yhat, lambda = lambda)
    class(result) <- "dyntermstrc_yields"
    result
  }



estimateyieldcurve <- function(method, y, m, beta, lambda, objfct, grad_objfct, constraints, constrOptimOptions) {
  if(method=="dl") {
    opt_result <- constrOptim(theta = beta,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              lambda, m, y)
  } else {
    opt_result <- constrOptim(theta = beta,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              m,y)
  }     
}

findstartparamyields <- function(y,m, method, tauconstr, control = list(), outer.iterations = 30, outer.eps = 1e-04)
  {
    epsConst <- 0.0001 ## ensures that starting value for constrOptim can not be on the boundary

    if(method=="dl"){
      startparam = rep(1,3)
      tau = NULL
      fmin = NULL
      optind = NULL
    }
    
    if(method=="ns"){
      tau <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
      fmin <- rep(NA, length(tau))
      lsbeta <- matrix(nrow = length(tau), ncol = 4)

      ui <- rbind(c(1,0,0),                 # beta0 > 0
                  c(1,1,0))                 # beta0 + beta1 > 0
      ci <- c(0,0)
      
      for (i in 1:length(tau)){
        X <- cbind(rep(1,length(y)),
                   ((1 - exp(-m/tau[i]))/(m/tau[i])),
                   (((1 - exp(-m/tau[i]))/(m/tau[i])) - exp(-m/tau[i])))

        lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
        beta <- c(lsparam[1:3],tau[i])

        ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0)
        if(beta[1]>0 && ((beta[1]+beta[2])>0)){
          fmin[i] <- objfct_ns(beta, m, y)
        } else {
          ## switch to constrOptim if OLS violates constraints
          lsparam <- constrOptim(theta = rep(1,3), # start parameters for D/L, objective function is convex
                                 f = objfct_ns_grid,
                                 grad = grad_ns_grid,
                                 ui = ui,
                                 ci = ci,
                                 mu = 1e-04,
                                 control = control,
                                 method = "BFGS",
                                 outer.iterations = outer.iterations,
                                 outer.eps = outer.eps,
                                 tau[i], m, y) ## additional inputs for f and grad
          
          beta <- c(lsparam$par,tau[i])
          fmin[i] <- objfct_ns(beta, m, y)
        }
        lsbeta[i,] <- beta 
      }
      optind <- which(fmin == min(fmin))
      startparam <- lsbeta[optind,]
    }

    
    if(method=="sv"){
      tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
      tau2 <- tau1
      tau <- cbind(tau1, tau2)
      fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
      ## fminsolver <- matrix(nrow = length(tau1), ncol = length(tau2)) ## DEBUG
      lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
      
      ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                  c(1,1,0,0))                 # beta0 + beta1 > 0
      ci <- c(0,0)
      
      for (i in 1:length(tau1))
        {
          for (j in 1:length(tau2))
            {
              ## minimum tau distance constraint
              if(tau1[i] + tauconstr[4] < tau2[j]) {
                
                ## reparametrize to avoid nonsingular matrix
                if (i == j){
                  X <- cbind(rep(1,length(y)),
                             ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                             - exp(-m/tau1[i]))
                  
                  lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
                  beta <- c(lsparam[1],lsparam[2]-lsparam[3],lsparam[3]/2,tau1[i],lsparam[3]/2,tau2[j])
                } else
                {
                  X <- cbind(rep(1,length(y)),
                             ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                             (((1 - exp(-m/tau1[i]))/(m/tau1[i])) - exp(-m/tau1[i])),
                             (((1 - exp(-m/tau2[j]))/(m/tau2[j])) - exp(-m/tau2[j])))
                  
                  lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
                  beta <- c(lsparam[1:3],tau1[i],lsparam[4],tau2[j])
                }
                ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0, tau distance)
                ##if(beta[1]>0 && ((beta[1]+beta[2])>0  ) && (-tau1[i] + tau2[j]) > tauconstr[4]){
                if(beta[1]>0 && ((beta[1]+beta[2])>0)){
                  fmin[i,j] <- objfct_sv(beta, m, y)
                } else {
                  ## switch to constrOptim if OLS violates constraints
                  lsparam <- constrOptim(theta = rep(0.01,4),
                                         f = objfct_sv_grid,
                                         grad = grad_sv_grid,
                                         ui = ui,
                                         ci = ci,
                                         mu = 1e-04,
                                         control = control,
                                         method = "BFGS",
                                         outer.iterations = outer.iterations,
                                         outer.eps = outer.eps,
                                         c(tau1[i], tau2[j]), m, y) ## additional inputs for f and grad
                  
                  beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
                  fmin[i,j] <- lsparam$value
                }
                lsbeta[(i-1)*length(tau1)+j,] <- beta
              }
            } 
        }
      optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
      startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]
    }
    
    if(method=="asv"){
      tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
      tau2 <- tau1
      tau <- cbind(tau1, tau2)
      fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
      lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
      
      ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                  c(1,1,0,0))                 # beta0 + beta1 > 0
      ci <- c(0,0)
      
      for (i in 1:length(tau1))
        {
          for (j in 1:length(tau2))
            {
              if(tau1[i] < tau2[j]) {
                
                X <- cbind(rep(1,length(y)),
                           ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                           (((1 - exp(-m/tau1[i]))/(m/tau1[i])) - exp(-m/tau1[i])),
                           (((1 - exp(-m/tau2[j]))/(m/tau2[j])) - exp(-2*m/tau2[j])))
                
                lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
                beta <- c(lsparam[1:3],tau1[i],lsparam[4],tau2[j])
                
                ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0, tau distance)
                ##if(beta[1]>0 && ((beta[1]+beta[2])>0) && (-tau1[i] + tau2[j]) > 0){
                if(beta[1]>0 && ((beta[1]+beta[2])>0)){
                  fmin[i,j] <- objfct_sv(beta, m, y)
                } else {
                  ##print(paste("OLS violated constraints, solving with constrOptim, tau_1 =",tau1[i],"and tau_2 = ", tau2[j])) ## DEBUG
                  ## switch to constrOptim if OLS violates constraints
                  lsparam <- constrOptim(theta = rep(0.01,4),
                                         f = objfct_asv_grid,
                                         grad = grad_asv_grid,
                                         ui = ui,
                                         ci = ci,
                                         mu = 1e-04,
                                         control = control,
                                         method = "BFGS",
                                         outer.iterations = outer.iterations,
                                         outer.eps = outer.eps,
                                         c(tau1[i], tau2[j]), m, y) ## additional inputs for f and grad
                  
                  beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
                  fmin[i,j] <- lsparam$value
                }
                lsbeta[(i-1)*length(tau1)+j,] <- beta
              }
            } 
        }
      optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
      startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]
    }
    
    result <- list(startparam = startparam, tau = tau, fmin = fmin, optind = optind)
    class(result) <- "spsearch"
    result
  }

