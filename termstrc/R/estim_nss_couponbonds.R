##########################################################################
### Nelson/Siegel-type yield curve estimation method for 'couponbonds' ###
##########################################################################

estim_nss <- function(dataset, ...) UseMethod("estim_nss")

estim_nss.couponbonds <- function(dataset,                  # dataset (static)
                                  group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
                                  matrange="all",            # maturity range in years c(min, max) 
                                  method="ns",
                                  startparam=NULL,           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
                                                             # otherwise globally optimal parameters are searched automatically
                                  lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                  tauconstr = NULL,          # constraints for tau parameter grid
                                  constrOptimOptions = list(control = list(maxit = 2000), outer.iterations = 200, outer.eps = 1e-04),...
           ) {

  ## data preprocessing
  prepro <- prepro_bond(group=group,bonddata=dataset,matrange=matrange)

  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
  
  ## automatically determine globally optimal start parameters
  spsearch <- list()
  length(spsearch) <- n_group

   ## default tau constraints (if not specified by user)
  if (is.null(tauconstr)) {
    tauconstr <- list()
    length(tauconstr) <- n_group
    names(tauconstr) <- group
    for (k in sgroup){
        tauconstr[[k]] <- c(min(m[[k]][1,]), max(m[[k]]), 0.2, 0.5)
        names(tauconstr[[k]]) <- c("tau_min", "tau_max", "gridstepsize", "tau_distance")
        if (method == "asv") {tauconstr[[k]][4] = 0}
        if (method == "ns") {tauconstr[[k]] = tauconstr[[k]][1:3]}
    }
    if (method!="dl"){
      print("The following constraints are used for the tau parameters:")
      print(tauconstr)}
  }
  
  if(is.null(startparam)){
    startparam <- matrix(ncol = 6, nrow = n_group)
    
    colnames(startparam) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
    
    if (method == "dl") {startparam <- startparam[,1:3, drop=FALSE]}
    if (method == "ns") {startparam <- startparam[,1:4, drop=FALSE]}

    for (k in sgroup){
      ## check if grid size is too large
      if((tauconstr[[k]][1] + tauconstr[[k]][3]) > (tauconstr[[k]][2] - tauconstr[[k]][3])) warning("Grid step size is too large!")
      if((tauconstr[[k]][1] + 3*tauconstr[[k]][3]) > tauconstr[[k]][2]) warning("Grid step size is too large!")
         
      print(paste("Searching startparameters for ", group[k]))
      spsearch[[k]] <- findstartparambonds(p[[k]],m[[k]],cf[[k]], duration[[k]][,3],
                                            method, tauconstr[[k]])
      startparam[k,] <- spsearch[[k]]$startparam 
      print(startparam[k,])
    }
  }
  
  rownames(startparam) <- group
  
  ## objective function (weighted price error minimization) 
  objfct <- get_objfct_bonds(method)

  ## gradient objective function
  grad_objfct <- get_grad_objfct_bonds(method)
  
  ## calculate optimal parameter vectors
  constraints <- list()
  for (k in sgroup){
    constraints[[k]] <- get_constraints(method, tauconstr[[k]])
  }
  opt_result <- list()
  for (k in sgroup){
    opt_result[[k]] <- estimatezcyieldcurve(method, startparam[k,], lambda, objfct, grad_objfct, constraints[[k]],
                                            constrOptimOptions, m[[k]], cf[[k]], duration[[k]][,3], p[[k]]) 
  }

  ## data post processing 
  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,method,lambda)
  
  ## return list of results 
  result <- list(group=group,                   # e.g. countries, rating classes
                 matrange=matrange,             # maturity range of bonds
                 method=method,                 # estimation method
                 startparam=startparam,         # calculated startparameters
                 n_group=n_group,               # number of groups,
                 lambda=lambda,                 # lambda parameter for dl
                 spsearch = spsearch,           # detailed data from start param search
                 spot=postpro$zcy_curves,       # zero coupon yield curves
                 spread=postpro$s_curves,       # spread curves
                 forward=postpro$fwr_curves,    # forward rate curves
                 discount=postpro$df_curves,    # discount factor curves
                 expoints=postpro$expoints,     # extrapolation points
       		 cf=cf,                         # cashflow matrix
                 m=m,                           # maturity matrix
                 duration=duration,             # duration, modified duration, weights
                 p=p,                           # dirty prices        
                 phat=postpro$phat,             # estimated dirty prices         
                 perrors=postpro$perrors,       # price errors
                 ac=ac,                         # accrued interest
                 y=y,                           # maturities and yields
                 yhat=postpro$yhat,             # estimated yields
                 yerrors=postpro$yerrors,       # yield errors
                 opt_result=opt_result          # optimisation results           
                 )
              
  for ( i in 7:length(result)) names(result[[i]]) <- group
  class(result) <- "termstrc_nss"
  result
}

### Estimate zero-coupon yield curve

estimatezcyieldcurve <- function(method, startparam, lambda, objfct, grad_objfct, constraints, constrOptimOptions, m, cf, weights, p) {

  if(method=="dl") {
      opt_result <- constrOptim(theta = startparam,
                                f = objfct,
                                grad = grad_objfct,
                                ui = constraints$ui,
                                ci = constraints$ci,
                                mu = 1e-04,
                                control = constrOptimOptions$control,
                                method = "BFGS",
                                outer.iterations = constrOptimOptions$outer.iterations,
                                outer.eps = constrOptimOptions$outer.eps,
                                lambda, m, cf, weights, p)
    } else {
      opt_result <- constrOptim(theta = startparam,
                                f = objfct,
                                grad = grad_objfct,
                                ui = constraints$ui,
                                ci = constraints$ci,
                                mu = 1e-04,
                                control = constrOptimOptions$control,
                                method = "BFGS",
                                outer.iterations = constrOptimOptions$outer.iterations,
                                outer.eps = constrOptimOptions$outer.eps,
                                m, cf, weights, p)
    }
    opt_result
}

### Start parameter search routine for bond data

findstartparambonds <- function(p,m,cf, weights, method, tauconstr,
                                control = list(), outer.iterations = 30, outer.eps = 1e-04) {

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
      
      lsparam <- constrOptim(theta = rep(1,3), # start parameters for D/L, objective function is convex
                               f = objfct_ns_bonds_grid,
                               grad = grad_ns_bonds_grid,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = control,
                               method = "BFGS",
                               outer.iterations = outer.iterations,
                               outer.eps = outer.eps,
                               tau[i], m, cf, weights, p) ## additional inputs for f and grad
      beta <- c(lsparam$par,tau[i])
      fmin[i] <- lsparam$value
      lsbeta[i,] <- beta 
    }
    optind <- which(fmin == min(fmin, na.rm = TRUE))
    startparam <- lsbeta[optind,]
  }
    
  if(method=="sv") {

    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)
      
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    for (i in 1:length(tau1))
      {
        #print(i) # DEBUG
        for (j in 1:length(tau2))
          {
            
            if(tau1[i] + tauconstr[4] < tau2[j]) {
              #print(j) # DEBUG
              
              lsparam <- constrOptim(theta = rep(0.01,4),
                                     f = objfct_sv_bonds_grid,
                                     grad = grad_sv_bonds_grid,
                                     ui = ui,
                                     ci = ci,
                                     mu = 1e-04,
                                     control = control,
                                     method = "BFGS",
                                     outer.iterations = outer.iterations,
                                     outer.eps = outer.eps,
                                     c(tau1[i], tau2[j]), m, cf, weights, p) ## additional inputs for f and grad
              
              beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
              
              fmin[i,j] <- lsparam$value
              lsbeta[(i-1)*length(tau1)+j,] <- beta
            }
          }
      }
    
    optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
    startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]    
  }

  if(method=="asv") {
    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)
      
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    for (i in 1:length(tau1))
      {
        for (j in 1:length(tau2))
          {
            
            if(tau1[i] < tau2[j]) {
              
              lsparam <- constrOptim(theta = rep(0.01,4),
                                     f = objfct_asv_bonds_grid,
                                     grad = grad_asv_bonds_grid,
                                     ui = ui,
                                     ci = ci,
                                     mu = 1e-04,
                                     control = control,
                                     method = "BFGS",
                                     outer.iterations = outer.iterations,
                                     outer.eps = outer.eps,
                                     c(tau1[i], tau2[j]), m, cf, weights, p) ## additional inputs for f and grad
              
              beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
              
              fmin[i,j] <- lsparam$value
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

### Startparameter grid search plots

plot.spsearch <- function(x, main = "Start parameter search", rgl = TRUE, ...) {

  if(is.matrix(x$tau)){
    image(x$tau[,1],x$tau[,2],log(x$fmin),xlab = expression(tau[1]), ylab = expression(tau[2]),main = "Log(Objective function)")
      contour(x$tau[,1],x$tau[,2],log(x$fmin),nlevels=10,xlab = expression(tau[1]), ylab = expression(tau[2]),main = "Log(Objective function)", add = TRUE)
      points(x$tau[x$optind[1],1],x$tau[x$optind[2],2],pch = 10, col = "steelblue")
      if (rgl) {
        open3d()
        persp3d(x$tau[,1], x$tau[,2], log(x$fmin), col = "green3", box = FALSE,xlab = "tau_1", ylab = "tau_2", zlab = "Log(Objective function)")
        points3d(x$tau[x$optind[1],1],x$tau[x$optind[2],2],min(log(x$fmin), na.rm = TRUE), col = "red")
      }
      else {
        par(ask = TRUE)
        persp(x$tau[,1], x$tau[,2], log(x$fmin), col = "green3", box = TRUE, xlab = "tau_1", ylab = "tau_2", zlab = "Log(Objective function)",
              shade = TRUE, ticktype = "detailed", border = NA, cex.lab = 1, cex.axis = 0.7,  theta = 0, phi = 25, r = sqrt(3),
              d = 1, scale = TRUE, expand = 1, ltheta = 135, lphi = 0)
      }
  } else {
      plot(x$tau,log(x$fmin),xlab = expression(tau[1]), ylab = "Log(Objective function)", type = "l", main = main)
      points(x$tau[x$optind],log(x$fmin[x$optind]),pch = 10, col = "red")
  }
}
