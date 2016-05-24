#######################################################
##     Bayesian compound Poisson Linear Models       ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

bcplm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter / 2),
                   n.thin = max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims = 1000, n.report = 2, prior.beta.mean = NULL, 
                   prior.beta.var = NULL, bound.phi = 100, bound.p = c(1.01, 1.99), 
                   tune.iter = 5000, n.tune = floor(tune.iter/100),
                   basisGenerators = c("tp", "bsp", "sp2d"), ...) {

  call <- expand.call(match.call())  
  if (missing(data)) 
    data <- environment(formula)   
  link.power <- make.link.power(link)
  # identify smooth terms 
  tf <- terms.formula(formula, specials = 
                eval(call$basisGenerators, parent.frame(2)))
  n.f <- length(unlist(attr(tf, "specials")))
  is.cpglmm <- 0
  
  # this is for cpglmm 
  if (length(findbars(formula)) || n.f){
    is.cpglmm <- 1
    # create model frame and get factor list
    if (n.f) {
      call2 <- as.list(call)[-1]
      m <- match(c("formula", "data", "weights", "offset", 
                   "contrasts", "basisGenerators"), names(call2),0L)
      call2 <- call2[m]
      setup <- do.call(frFL, as.list(call2))
      fr <- setup$m$fr 
      FL <- setup$m$FL
    } else {  
      fr <- lmerFrames(call, formula, contrasts)
      FL <- lmerFactorList(formula, fr, 0L, 0L)      
    } 
    dm <- mkZt(FL, NULL)
  } else {   
    # this is for glm
    fr <- cpglm.mf(call, contrasts)    
  }
  # dimensions 
  n.obs <- NROW(fr$X)
  n.beta <- NCOL(fr$X)
  
  # default weights and offsets if NULL (for cpglmm)   
  if (is.null(fr$wts) || length(fr$wts) == 0) fr$wts <- as.double(rep(1, n.obs)) 
  if (is.null(fr$off) || length(fr$off) == 0) fr$off <- as.double(rep(0, n.obs))  
  
  # check arguments
  check.args.bcplm(call, n.beta, n.chains)
  
  # default prior mean and var if missing           
  if (is.null(prior.beta.mean)) prior.beta.mean <- rep(0, n.beta)			
  if (is.null(prior.beta.var)) prior.beta.var <- rep(10000, n.beta)
  
  # the dims slot in bcplm_input
  n.keep <- floor((n.iter - n.burnin)/n.thin)
  n.sims <- n.chains * n.keep
  n.report <- ifelse(n.report, floor(n.iter/n.report), 0)
  dims <- as.integer(c(n.obs, n.beta, sum(fr$Y > 0), 0, 0, n.beta + 2,
            n.chains, n.iter, n.burnin, n.thin, n.keep,
            n.sims, n.report, tune.iter, n.tune, n.beta + 2))  
  names(dims) <- c("n.obs", "n.beta", "n.pos", "n.term", "n.u", 
            "n.all", "n.chains", "n.iter", "n.burnin", "n.thin", 
            "n.keep", "n.sims", "n.report", "tune.iter", "n.tune", "n.mh")  
  if (is.cpglmm){
    ncol <- as.integer(unlist(lapply(dm$ST, ncol)))
    dims[c("n.term", "n.u")] <- unname(dm$dd[c("nt", "q")])
    dims["n.all"] <- as.integer(dims["n.all"] + dims["n.u"] + sum(ncol^2))
    dims["n.mh"] <- as.integer(n.beta + dims["n.u"] + 2)
  }
  
  # proposal sd's
  mh.sd <- rep(1, dims["n.mh"])               
  
  # generate initial values if necessary
  if (is.null(inits)) {
    if (!is.cpglmm) dm <- NULL
    tmp <- bcplm.init(fr, link.power, n.chains, bound.p, dm)
    inits <- tmp$inits  
    # update proposal covariance matrix
    mh.sd[1:n.beta] <- sqrt(diag(tmp$vbeta))
  } else{
    # check initial values
    check.inits.bcplm(inits, n.beta, dims["n.term"], n.chains, is.cpglmm)  
    inits <- lapply(inits, function(x) 
                c(x$beta, x$phi, x$p, x$u, unlist(lapply(x$Sigma, as.numeric))))
  }
  
  # input for the C functions  
  input <- new("bcplm_input", X = fr$X, y = as.double(fr$Y), 
             Zt = if (is.cpglmm) dm$Zt else as(matrix(0), "dgCMatrix"),
             ygt0 = as.integer(which(fr$Y > 0L) - 1),
             offset = as.double(fr$off), pWt= as.double(fr$wts),
             mu = double(n.obs), eta = double(n.obs),
             Xb = double(n.obs), Zu = double(n.obs),   
             inits = inits, fixef = as.double(inits[[1]][1:n.beta]),
             phi = as.double(inits[[1]][n.beta + 1]),
             p = as.double(inits[[1]][n.beta + 2]),
             u = as.double(inits[[1]][(n.beta + 3):(n.beta + 2 + dims["n.u"])]), 
             Sigma = if (!is.cpglmm) list() else lapply(dm$ST, function(x) x %*% t(x)),
             link.power = as.double(link.power),
             pbeta.mean = as.double(prior.beta.mean),
             pbeta.var = as.double(prior.beta.var),
             bound.phi = as.double(bound.phi),
             bound.p = as.double(bound.p),    
             mh.sd = as.double(mh.sd), dims = dims,  
             k = as.integer(0), cllik = double(1),
             Gp = if (!is.cpglmm) as.integer(0) else unname(dm$Gp),   
             ncol = if (!is.cpglmm) as.integer(0) else ncol,     
             nlev = if (!is.cpglmm) as.integer(0) else 
                    as.integer(sapply(FL$fl, function(x) length(levels(x)))),
             accept = double(n.beta + 2 + dims["n.u"]))
  
  # run MCMC
  sims.list <- .Call("bcplm_mcmc", input)  
  
  # set names
  xnm <- if (!is.cpglmm) dimnames(fr$X)[[2L]] else names(fr$fixef)
  parnm <- c(xnm, "phi", "p")
  if (is.cpglmm){
    unm <- paste("u", 1:dm$dd['q'], sep = "")
    snm <- sapply(1:length(dm$ST), function(x){
                    tm <- apply(expand.grid(1:ncol[x], 1:ncol[x]), 1,
                    paste, collapse = ",")
                    paste("Sigma", x, "[", tm, "]", sep = "")})
    parnm <- c(parnm, unm, snm)
  }
  sims.list <- lapply(sims.list, function(x){ 
                    dimnames(x) <- list(NULL, parnm)
                    return(x)})
  
  # coerce to mcmc.list                  
  sims.list <- as.mcmc.list(lapply(sims.list, as.mcmc))
  s <- summary(sims.list)
  # output simulation results 
  psd <- input@mh.sd
  psd <- list(beta = psd[1:n.beta], phi = psd[n.beta + 1],
              p = psd[n.beta + 2],
              u = if(is.cpglmm) psd[-(1:(n.beta + 2))] else list())
  
  # get estimated variance components
  if (is.cpglmm) {
    Sigma <- getSigmaList(dims["n.term"], s)
    Sigma <- lapply(1:length(Sigma), function(x){
                tmp <- Sigma[[x]]
                dimnames(tmp) <- dimnames(dm$ST[[x]])
                tmp
              })
  }
  ans <- new("bcplm", 
             dims = dims, sims.list = sims.list,
             link.power = link.power, call = call,
             formula = formula, model.frame = fr$mf,
             contrasts = contrasts, inits = inits,
             Zt = if(is.cpglmm) dm$Zt else as(matrix(0), "dgCMatrix"), 
             flist = if(is.cpglmm) dm$flist else list(),
             summary = s,
             Sigma = if(is.cpglmm) Sigma else list(),  
             prop.sd = psd)  
  return(ans)
}


# get summary of Sigma in a list format
getSigmaList <- function(nt, s){
  rn <- rownames(s[[1]])
  Sigma <- lapply(1:nt, function(xx){
    tmp <- s[[2]][grep(paste("^Sigma", xx, "\\[", sep = ""), rn), 3]
    matrix(tmp, sqrt(length(tmp)))
  })
  Sigma
}






