##########################################################################
## samples from the posterior distribution of a dynamic 1d IRT model
## a la Martin and Quinn (2002)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Assumes a local level model for the evolution of theta
##
##  y_{jkt}^* = -alpha_k + beta_k * theta_{jt} + epsilon_{jkt}
##  theta_{jt} ~ N(theta_{j(t-1)}, tau^2)
##
## Kevin Quinn
## 1/28/2008
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMCdynamicIRT1d" <- function(datamatrix, item.time.map,
                               theta.constraints=list(),
                               burnin=1000, mcmc=20000, thin=1,
                               verbose=0, seed=NA, theta.start=NA,
                               alpha.start=NA, beta.start=NA,
                               tau2.start=1,
                               a0=0, A0=.1, b0=0, B0=.1,
                               c0=-1, d0=-1, e0=0, E0=1,
                               store.ability=TRUE,
                               store.item=TRUE, ... 
                               ){

  
  datamatrix <- as.matrix(datamatrix)
  
  nitems <- ncol(datamatrix)
  nsubj  <- nrow(datamatrix)
  ntime <- max(item.time.map)
  
  ## checks
  check.offset(list(...))
  check.mcmc.parameters(burnin, mcmc, thin)

  
  if (nitems != length(item.time.map)){
    cat("Number of rows of datamatrix not equal to length of item.time.map\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)    
  }
  if (min(item.time.map) != 1){
    cat("Minimum value in item.time.map not equal to 1\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)    
  }
  if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) !=
     (nitems * nsubj)) {
    cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
    stop("Please check data and try MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }  
  
  
  if (A0 < 0){
    cat("Error: A0 (prior precision for alpha) is less than 0).\n")
    stop("Please respecify and try MCMCdynamicIRT1d() again.\n")
  }
  if (B0 < 0){
    cat("Error: B0 (prior precision for beta) is less than 0).\n")
    stop("Please respecify and try MCMCdynamicIRT1d() again.\n")
  }
  




  
  ## setup constraints on theta
  if(length(theta.constraints) != 0) {
    for (i in 1:length(theta.constraints)){
      theta.constraints[[i]] <-
        list(as.integer(1), theta.constraints[[i]][1])
    }
  }
  holder <- build.factor.constraints(theta.constraints, t(datamatrix), nsubj, 1)
  theta.eq.constraints <- holder[[1]]
  theta.ineq.constraints <- holder[[2]]
  ##subject.names <- holder[[3]]
  ## names
  item.names <- colnames(datamatrix)
  if (is.null(item.names)){
    item.names <- paste("item", 1:nitems, sep="")
  }
  




  
  ## starting values for theta error checking
  theta.start <- factor.score.start.check(theta.start, datamatrix,
                                          matrix(0,1,1), matrix(1,1,1),
                                          theta.eq.constraints,
                                          theta.ineq.constraints, 1)
  
  ## starting values for (alpha, beta)
  ab.starts <- matrix(NA, nitems, 2)
  for (i in 1:nitems){
    local.y <-  datamatrix[,i]
    ##local.y[local.y==9] <- NA
    if (length(na.omit(local.y)) <= 1){
      ab.starts[i,] <- c(0, 10)
    }
    else if (var(na.omit(local.y))==0){
      ab.starts[i,] <- c(0,10)
      
    }
    else {
      ab.starts[i,] <- coef(suppressWarnings(glm(local.y~theta.start,
                                                 family=binomial(link="probit"),
                                                 control=glm.control(
                                                   maxit=8, epsilon=1e-3)
                                                 )))
    }
  }
  ab.starts[,1] <- -1 * ab.starts[,1] # make this into a difficulty param
  
  ## starting values for alpha and beta error checking
  if (is.na(alpha.start)) {
    alpha.start <- ab.starts[,1]
  }
  else if(is.null(dim(alpha.start))) {
    alpha.start <- alpha.start * matrix(1,nitems,1)  
  }
  else if((dim(alpha.start)[1] != nitems) || (dim(alpha.start)[2] != 1)) {
    cat("Error: Starting value for alpha not conformable.\n")
    stop("Please respecify and call MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }      
  if (is.na(beta.start)) {
    beta.start <- ab.starts[,2]
  }
  else if(is.null(dim(beta.start))) {
    beta.start <- beta.start * matrix(1,nitems,1)  
  }
  else if((dim(beta.start)[1] != nitems) || (dim(beta.start)[2] != 1)) {
    cat("Error: Starting value for beta not conformable.\n")
    stop("Please respecify and call MCMCdynamicIRT1d() again.\n",
         call.=FALSE)
  }    
  



  ## generate time-specific theta information and create extended theta.start
  subject.names <- NULL
  theta.start.new <- NULL
  ## theta.info has:
  ## col1: subj ID, col2: #time periods, col3: offset (first term C indexing)
  theta.info <- matrix(NA, nsubj, 3)
  for (s in 1:nsubj){
    namestub <- rownames(datamatrix)[s]
    theta.info[s,1] <- s
    count <- 0
    holder <- NULL
    for (i in 1:nitems){
      if (!is.na(datamatrix[s,i])){
        holder <- c(holder, item.time.map[i])
      }
    }
    minholder <- min(holder)
    maxholder <- max(holder)    
    theta.info[s,2] <- maxholder - minholder + 1
    theta.info[s,3] <- minholder - 1
    theta.start.new <- c(theta.start.new, rep(theta.start[s], theta.info[s,2]))
    subject.names <- c(subject.names,
                       paste(namestub, ".t", minholder:maxholder, sep=""))
  }
  nthetas <- length(subject.names)
  theta.start <- theta.start.new


  
  if (length(c0) < nsubj){
    c0 <- array(c0, nsubj)
  }
  if (length(d0) < nsubj){
    d0 <- array(d0, nsubj)
  }
  if (length(tau2.start) < nsubj){
    tau2.start <- array(tau2.start, nsubj)
  }
  if (length(e0) < nsubj){
    e0 <- array(e0, nsubj)
  }
  if (length(E0) < nsubj){
    E0 <- array(E0, nsubj)
  }
  E0inv <- 1/E0

  
  subject.names.short <- rownames(datamatrix)

  
  ## convert datamatrix into a sparse format where the missing values are
  ## not explicitly represented
  ##
  ## 1st column: subject index (C indexing)
  ## 2nd column: item index (C indexing)
  ## 3rd column: vote 
  data.sparse.si <- NULL
  for (i in 1:nsubj){
    for (j in 1:nitems){
      if (!is.na(datamatrix[i,j])){
        data.sparse.si <- rbind(data.sparse.si, c(i-1, j-1, datamatrix[i,j]))
      }
    }
  }
  ## 1st column: item index (C indexing)
  ## 2nd column: subject index (C indexing)
  ## 3rd column: vote 
##  data.sparse.is <- NULL
##  for (i in 1:nitems){
##    for (j in 1:nsubj){
##      if (!is.na(datamatrix[j,i])){
##        data.sparse.is <- rbind(data.sparse.is, c(i-1, j-1, datamatrix[j,i]))
##      }
##    }
##  }
  
  rm(datamatrix)




  

  ## define storage matrix for posterior theta draws
  thetadraws <- matrix(0, mcmc/thin, nthetas)
  if (store.ability != TRUE){
    thetadraws <- matrix(0, 1, 1)
  }
  alphadraws <- matrix(1, mcmc/thin, nitems)
  betadraws  <- matrix(2, mcmc/thin, nitems)
  if (store.item != TRUE){
    alphadraws <- matrix(1, 1, 1)
    betadraws  <- matrix(2, 1, 1)
  }
  tau2draws <- matrix(0, mcmc/thin, nsubj)


  ## seeds
  seeds <- form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
 ## print(theta.eq.constraints)
 ## print(theta.ineq.constraints)


#  return(list(theta=theta.start, alpha=alpha.start, beta=beta.start))

  
  ## call C++ code to draw sample
  posterior <- .C("MCMCdynamicIRT1d",
                  
                  thetadata = as.double(thetadraws),
                  thetarow = as.integer(nrow(thetadraws)),
                  thetacol = as.integer(ncol(thetadraws)),
                  
                  alphadata = as.double(alphadraws),
                  alpharow = as.integer(nrow(alphadraws)),
                  alphacol = as.integer(ncol(alphadraws)),                  

                  betadata = as.double(betadraws),
                  betarow = as.integer(nrow(betadraws)),
                  betacol = as.integer(ncol(betadraws)),

                  tau2data = as.double(tau2draws),
                  tau2row = as.integer(nrow(tau2draws)),
                  tau2col = as.integer(ncol(tau2draws)),
                  
                  nsubj = as.integer(nsubj),
                  nitems = as.integer(nitems),
                  ntime = as.integer(ntime),
                  
                  Ydata = as.integer(data.sparse.si),
                  Yrow = as.integer(nrow(data.sparse.si)),
                  Ycol = as.integer(ncol(data.sparse.si)),

                  IT = as.integer(item.time.map-1),
                  ITlength = as.integer(length(item.time.map)),
                  
                  burnin = as.integer(burnin),
                  mcmc = as.integer(mcmc),
                  thin = as.integer(thin),
                  
                  lecuyer = as.integer(lecuyer),
                  seedarray = as.integer(seed.array),
                  lecuyerstream = as.integer(lecuyer.stream),
                  
                  verbose = as.integer(verbose),

                  thetastartdata = as.double(theta.start),
                  thetastartlength = as.integer(length(theta.start)),

                  thetainfo = as.integer(theta.info),
                  thetainforow = as.integer(nrow(theta.info)),
                  thetainfocol = as.integer(ncol(theta.info)),
                  
                  astartdata = as.double(alpha.start),
                  astartlength = as.integer(length(alpha.start)),
 
                  bstartdata = as.double(beta.start),
                  bstartlength = as.integer(length(beta.start)),

                  tau2data = as.double(tau2.start),
                  tau2length = as.integer(length(tau2.start)),

                  c0data = as.double(c0),
                  c0length = as.integer(length(c0)),

                  d0data = as.double(d0),
                  d0length = as.integer(length(d0)),
                  
                  a0data = as.double(a0),
                  A0data = as.double(A0),
                  
                  b0data = as.double(b0),
                  B0data = as.double(B0),

                  e0data = as.double(e0),
                  E0invdata = as.double(E0inv),

                  thetaeqdata = as.double(theta.eq.constraints),
                  thetaeqrow = as.integer(nrow(theta.eq.constraints)),
                  thetaeqcol = as.integer(ncol(theta.eq.constraints)),

                  thetaineqdata = as.double(theta.ineq.constraints),
                  thetaineqrow = as.integer(nrow(theta.ineq.constraints)),
                  thetaineqcol = as.integer(ncol(theta.ineq.constraints)),

                  storei = as.integer(store.item),
                  storea = as.integer(store.ability),
                  PACKAGE="MCMCpack"
                  )
  
  



  tau2samp <- matrix(posterior$tau2data,
                      posterior$tau2row,
                      posterior$tau2col)
  colnames(tau2samp) <- paste("tau2.", subject.names.short, sep="")



  if (store.item & store.ability){
    thetasamp <- matrix(posterior$thetadata,
                        posterior$thetarow,
                        posterior$thetacol) 
    colnames(thetasamp) <- paste("theta.", subject.names, sep="")
    
    alphasamp <- matrix(posterior$alphadata,
                        posterior$alpharow,
                        posterior$alphacol)
    colnames(alphasamp) <- paste("alpha.", item.names, sep="")
    
    
    betasamp <- matrix(posterior$betadata,
                       posterior$betarow,
                       posterior$betacol)
    colnames(betasamp) <- paste("beta.", item.names, sep="")
    
    outmat <- mcmc(cbind(thetasamp, alphasamp, betasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else if (store.item){
    alphasamp <- matrix(posterior$alphadata,
                        posterior$alpharow,
                        posterior$alphacol)
    colnames(alphasamp) <- paste("alpha.", item.names, sep="")
    
    
    betasamp <- matrix(posterior$betadata,
                       posterior$betarow,
                       posterior$betacol)
    colnames(betasamp) <- paste("beta.", item.names, sep="")
    
    outmat <- mcmc(cbind(alphasamp, betasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else if (store.ability){
    thetasamp <- matrix(posterior$thetadata,
                        posterior$thetarow,
                        posterior$thetacol) 
    colnames(thetasamp) <- paste("theta.", subject.names, sep="")

    outmat <- mcmc(cbind(thetasamp, tau2samp),
                   start=1, end=mcmc, thin=thin)
  }
  else {
    outmat <- mcmc(tau2samp,
                   start=1, end=mcmc, thin=thin)
  }
  
  return(outmat)

  
}


