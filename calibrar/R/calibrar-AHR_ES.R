
.continueEvolution = function(opt, control) {
  out = (opt$gen <= (control$maxgen - 1)) & (opt$step >= control$convergence)
#   out = (opt$gen <= control$maxit)
# reltol * (abs(val) + reltol)
  return(out)
}

# Initialize population ---------------------------------------------------

.calculateSigma = function(x, control) {
  sigma = control$sigma
  if(is.null(control$sigma)) sigma = rep(NA, length(x))
  out = x^2/12
  out[!is.finite(out)] = 1
  sigma[is.na(sigma)] = out[is.na(sigma)]
  return(sigma)
}

.calculateRange = function(lower, upper) {
  out = upper - lower
  out[!is.finite(out)] = 1
  return(out)
}

.newOpt = function(par, lower, upper, control) {
  
  opt = list()
  opt$gen       = 0
  opt$npar      = length(par)
  opt$nvar      = control$nvar
  opt$seed      = control$popsize
  opt$selection = control$selection 
  opt$range     = .calculateRange(lower, upper)
  opt$MU        = as.numeric(par)
  opt$SIGMA	    = .calculateSigma(opt$range, control)
  opt$lower     = lower
  opt$upper     = upper
  opt$ps		    = numeric(opt$npar)
  opt$pc		    = numeric(opt$npar)
  opt$ages	    = numeric(control$popsize)
  opt$step	    = control$step   
  opt$mu		    = matrix(opt$MU, nrow=opt$npar, ncol=opt$nvar)
  opt$sigma	    = matrix(0, nrow=opt$npar, ncol=opt$nvar)
  
  opt$SD        = opt$step*sqrt(opt$SIGMA)
  
  opt$parents   = ceiling(opt$selection*opt$seed)
  opt$w.rec     = .getRecombinationWeights(parents=opt$parents, method=control$method)
  opt$mu.eff    = 1/sum(opt$w.rec*opt$w.rec)
  opt$cc        = 4/(opt$npar+4)
  opt$chiN      = sqrt(opt$npar)*(1-1/(4*opt$npar)+1/(21*opt$npar^2))
  
  opt$cs        = (opt$mu.eff+2)/(opt$npar+opt$mu.eff+3)
  opt$D         = 1 + 2*max(0, sqrt((opt$mu.eff-1)/(opt$npar+1))-1) + opt$cs
  opt$mu.cov    = opt$mu.eff
  opt$c.cov     = (1/opt$mu.cov)*(2/((opt$npar + sqrt(2))^2)) + 
    (1-1/opt$mu.cov)*min(1,(2*opt$mu.eff-1)/((opt$npar+2)^2+opt$mu.eff))
  
  opt$aggFn     = match.fun(control$aggFn)
  opt$weights   = control$weights
  opt$useCV     = control$useCV
  opt$alpha     = control$alpha
  
  opt$control   = control
  
  class(opt) = c("restart", class(opt))
  return(opt)
  
}

.createPopulation = function(opt) {
  
  out = rtnorm2(n=opt$seed, mean=opt$MU, sd=opt$SD, lower=opt$lower, upper=opt$upper)
  return(out)

}

.getRecombinationWeights = function(parents, method) {
  out = log(parents+1)-log(1:parents)
  out = out/sum(out)
  return(out)
}


# Calculate Fitness -------------------------------------------------------


.calculateFitness = function(opt, fn) {
  
  fn = match.fun(fn)
  pop = opt$pop
  parallel = opt$control$parallel
  run      = opt$control$run
  
  path.tmp = getwd()               # get the current path
  on.exit(setwd(path.tmp))         # back to the original path after execution
  
  if(isTRUE(parallel)) {
    # optimize parallel execution, reduce data transfer
    
    FITNESS  =  foreach(i=0:(opt$seed-1), .combine=rbind, .verbose=FALSE, .inorder=FALSE) %dopar% {
      
      work.dir = .setWorkDir(run, i)    # set the 'individual' current directory
      
      Fitness = fn(pop[, i+1]) 
      Fitness = c(i+1, Fitness)
      Fitness
      
    }
    
    FITNESS = FITNESS[sort(FITNESS[,1], index.return=TRUE)$ix,][,-1, drop=FALSE]
    
  } else {
    
    FITNESS	=	NULL
    
    for(i in 0:(opt$seed-1)) {
      
      work.dir = .setWorkDir(run, i)    # set the 'individual' current directory

      Fitness = fn(pop[,i+1])
      FITNESS	=	rbind(FITNESS, Fitness)
      
    }
    
  }
  
  return(FITNESS)
}


# Selection ---------------------------------------------------------------


.selection = function(opt) {
  
  pop     = opt$pop
  fitness = opt$fitness
  
  fitness.global = .globalFitness(fitness=fitness, opt=opt)
  
  p = sort(fitness.global,index.return=TRUE)
  
  supsG			= p$ix[seq_len(opt$parents)]
  .best     = p$ix[1]
  
  best = list(BEST       = pop[, .best],
              fit	       = fitness[.best, ],
              fit.global = fitness.global[.best])
  
  supsL = apply(fitness[supsG, ,drop=FALSE], 2, FUN = function(x) sort(x, index.return=TRUE)$ix)
  
  return(list(supsG=supsG, supsL=supsL, best=best))
  
}


# Recombination -----------------------------------------------------------

.norma = function(x) {
  
  x[x==0] = 1E-20
  out = x/sum(x, na.rm=TRUE)
  return(out)
  
}

.w.oi=function(x, b=4) {
  
  n       =	ncol(x)
  if(n==1) return(matrix(1, nrow=nrow(x), ncol=1))
  
  cv.min  = 0.9*(apply(x, 2, min, na.rm=TRUE) + 1e-20)
  cv.max  = 1.1*(apply(x, 2, max, na.rm=TRUE) + 1e-20)
  out     = t(x)
  out     = ((cv.max - out)/(cv.max-cv.min))^b
  out     = t(out/rowSums(out))
  out	    = t(apply(out, 1, .norma))
  return(out)
  
}


.calculateOptimalWeights = function(opt) {
  
  pop   = opt$pop[, opt$selected$supsG]
  supsL = opt$selected$supsL
  
  w.rec = opt$w.rec
  alpha = opt$alpha
  
  opt.ind	= array(NA, dim=c(nrow(pop), ncol(supsL)))
  opt.sd	= array(NA, dim=c(nrow(pop), ncol(supsL)))
  
  for(i in seq_len(ncol(supsL))) {
    opt.ind[, i]       = apply(pop[, supsL[, i]], 1, weighted.mean, w=w.rec)
    opt.var	       = apply(pop[, supsL[, i]]^2, 1, weighted.mean, w=w.rec) - opt.ind[, i]^2
    opt.var[opt.var<0] = 0
    opt.sd[, i]        = sqrt(opt.var)
    
  }
  
  mu.new		= (1-alpha)*opt$mu + alpha*opt.ind
  s.new			= (1-alpha)*(opt$mu^2+opt$sigma^2) + alpha*(opt.ind^2+ opt.sd^2) - mu.new^2
  s.new[s.new<0]        = 0 # s.new is always positive, this is a correction for rounding errors 
  sigma.new             = sqrt(s.new)
  ww.new		= if(isTRUE(opt$useCV)) .w.oi(sigma.new/opt$range) else .w.oi(sigma.new)
  
  ww.rec = array(w.rec[supsL], dim=dim(supsL))
  W      = ww.new %*% t(ww.rec)

  opt$W     = W
  opt$w.new = ww.new
  opt$mu    = mu.new
  opt$sigma = sigma.new
  opt$best  = opt$selected$best # updateBestIndividual()
  
  opt$MU.eff  = 1/rowSums(W*W)
  opt$mu.eff  = max(mean(opt$MU.eff), 1) # mu.eff is always greater or equal to 1, avoiding rounding errors
    
  return(opt)
  
}


# .updatePopulation = function(opt) {
#   
#   opt = within(opt, {
#     
#     supsG       = selected$supsG
#     
#     cs          = (mu.eff+2)/(npar+mu.eff+3)
#     D           = 1 + 2*max(0, sqrt((mu.eff-1)/(npar+1))-1) + cs
#     mu.cov      = mu.eff
#     c.cov       = (1/mu.cov)*(2/((npar+sqrt(2))^2))+(1-1/mu.cov)*min(1,(2*mu.eff-1)/((npar+2)^2+mu.eff))
#     
#     MU.new      = rowSums(W * pop[, supsG, drop=FALSE])
#     
#     pc          = (1-cc)*pc + sqrt(cc*(2-cc))*sqrt(MU.eff)*(MU.new-MU)/step
#     ps          = (1-cs)*ps + sqrt(cs*(2-cs))*sqrt(MU.eff)*((MU.new-MU)/sqrt(SIGMA))/step
#     
#     SIGMA.sel   = rowSums(W * ((pop[, supsG, drop=FALSE] - MU)/step)^2)
#     SIGMA.new   = (1-c.cov)*SIGMA + (c.cov/mu.cov)*pc*pc + c.cov*(1-1/mu.cov)*SIGMA.sel
#     
#     step        = step*exp(cs*(sqrt(sum(ps*ps, na.rm=TRUE))/chiN-1)/D)
#   
#     MU          = MU.new
#     SIGMA       = SIGMA.new
#     
#     SD          = step*sqrt(SIGMA)
#     
#     rm(list=c("MU.new", "SIGMA.sel", "SIGMA.new"))
#     
#   })
#   
#   return(opt)
# 
# }

.updatePopulation = function(opt) {
  
    supsG       = opt$selected$supsG
    npar        = opt$npar
    MU          = opt$MU
    SIGMA       = opt$SIGMA
    step        = opt$step
    
    opt$cs      = (opt$mu.eff+2)/(npar+opt$mu.eff+3)
    opt$D       = 1 + 2*max(0, sqrt((opt$mu.eff-1)/(npar+1))-1) + opt$cs
    opt$mu.cov  = opt$mu.eff
    opt$c.cov   = (1/opt$mu.cov)*(2/((npar+sqrt(2))^2))+(1-1/opt$mu.cov)*
      min(1,(2*opt$mu.eff-1)/((npar+2)^2+opt$mu.eff))
    
    MU.new      = rowSums(opt$W * opt$pop[, supsG, drop=FALSE])
    
    opt$pc      = (1-opt$cc)*opt$pc + sqrt(opt$cc*(2-opt$cc))*sqrt(opt$MU.eff)*(MU.new-MU)/step
    opt$ps      = (1-opt$cs)*opt$ps + sqrt(opt$cs*(2-opt$cs))*sqrt(opt$MU.eff)*((MU.new-MU)/sqrt(SIGMA))/step
    
    SIGMA.sel   = rowSums(opt$W * ((opt$pop[, supsG, drop=FALSE] - MU)/step)^2)
    SIGMA.new   = (1-opt$c.cov)*SIGMA + (opt$c.cov/opt$mu.cov)*opt$pc*opt$pc + 
      opt$c.cov*(1-1/opt$mu.cov)*SIGMA.sel
    
    opt$step    = step*exp(opt$cs*(sqrt(sum(opt$ps*opt$ps, na.rm=TRUE))/opt$chiN-1)/opt$D)
    
    opt$MU      = MU.new
    opt$SIGMA   = SIGMA.new
    
    opt$SD      = opt$step*sqrt(opt$SIGMA)
    
  return(opt)
  
}


