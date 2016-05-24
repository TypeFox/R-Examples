##########################################################################
## sample from the posterior distribution of a one-dimensional item
## response theory model in R using linked C++ code in Scythe.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## ADM and KQ 1/23/2003
## updated extensively ADM & KQ 7/28/2004
## store.ability arg added KQ 1/27/2006
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCirt1d" <-
  function(datamatrix, theta.constraints=list(), burnin = 1000,
           mcmc = 20000, thin=1, verbose = 0, seed = NA,
           theta.start = NA, alpha.start = NA, beta.start = NA, t0 = 0,
           T0 = 1, ab0=0, AB0=.25, store.item = FALSE, store.ability=TRUE,
           drop.constant.items=TRUE, ... ) {
    
    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## check vote matrix and convert to work with C++ code
    if (drop.constant.items==TRUE){
      x.col.var <- apply(datamatrix, 2, var, na.rm=TRUE)
      keep.inds <- x.col.var>0
      keep.inds[is.na(keep.inds)] <- FALSE      
      datamatrix <- datamatrix[,keep.inds]
    }
    datamatrix <- as.matrix(datamatrix)   
    K <- ncol(datamatrix)   # cases, bills, items, etc
    J <- nrow(datamatrix)   # justices, legislators, subjects, etc
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (J*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and try MCMCirt1d() again.\n",
              call.=FALSE)
    }
    datamatrix[is.na(datamatrix)] <- 9   
    item.names <- colnames(datamatrix)
    subject.names <- rownames(datamatrix)
    
    ## setup constraints on theta
    if(length(theta.constraints) != 0) {
       for (i in 1:length(theta.constraints)){
          theta.constraints[[i]] <-
             list(as.integer(1), theta.constraints[[i]][1])
       }
    }
    holder <- build.factor.constraints(theta.constraints, t(datamatrix), J, 1)
    theta.eq.constraints <- holder[[1]]
    theta.ineq.constraints <- holder[[2]]
    subject.names <- holder[[3]]
    ## names
    item.names <- colnames(datamatrix)
    if (is.null(item.names)){
      item.names <- paste("item", 1:K, sep="")
    }

    ## prior for theta 
    holder <- form.mvn.prior(t0, T0, 1)
    t0 <- holder[[1]]
    T0 <- holder[[2]]

    ## prior for (alpha, beta)
    holder <- form.mvn.prior(ab0, AB0, 2)
    ab0 <- holder[[1]]
    AB0 <- holder[[2]]
    
    ## starting values for theta error checking
    theta.start <- factor.score.start.check(theta.start, datamatrix,
                                            t0, T0,
                                            theta.eq.constraints,
                                            theta.ineq.constraints, 1)
    
    ## starting values for (alpha, beta)
    ab.starts <- matrix(NA, K, 2)
    for (i in 1:K){
      local.y <-  datamatrix[,i]
      local.y[local.y==9] <- NA
      if (var(na.omit(local.y))==0){
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
      alpha.start <- alpha.start * matrix(1,K,1)  
    }
    else if((dim(alpha.start)[1] != K) || (dim(alpha.start)[2] != 1)) {
      cat("Error: Starting value for alpha not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
           call.=FALSE)
    }      
    if (is.na(beta.start)) {
      beta.start <- ab.starts[,2]
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)  
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Error: Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
              call.=FALSE)
    }    

    ## define holder for posterior sample
    if(store.item == FALSE & store.ability == TRUE) {
      sample <- matrix(data=0, mcmc/thin, J)
    }
    else if (store.item == TRUE & store.ability == FALSE){
      sample <- matrix(data=0, mcmc/thin, 2*K)
    }
    else if (store.item == TRUE & store.ability == TRUE){
      sample <- matrix(data=0, mcmc/thin, J + 2 * K)
    }
    else{
      cat("Error: store.item == FALSE & store.ability == FALSE.\n")
      stop("Please respecify and call MCMCirt1d() again.\n",
              call.=FALSE)      
    }

    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    # call C++ code to draw sample
    posterior <- .C("MCMCirt1d",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    Xdata = as.integer(datamatrix),
                    Xrow = as.integer(nrow(datamatrix)),
                    Xcol = as.integer(ncol(datamatrix)),    
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    thetastartdata = as.double(theta.start),
                    thetastartrow = as.integer(nrow(theta.start)),
                    thetastartcol = as.integer(ncol(theta.start)),
                    astartdata = as.double(alpha.start),
                    astartrow = as.integer(length(alpha.start)),
                    astartcol = as.integer(1),
                    bstartdata = as.double(beta.start),
                    bstartrow = as.integer(length(beta.start)),
                    bstartcol = as.integer(1),
                    t0 = as.double(t0),
                    T0 = as.double(T0),
                    ab0data = as.double(ab0),
                    ab0row = as.integer(nrow(ab0)),
                    ab0col = as.integer(ncol(ab0)),
                    AB0data = as.double(AB0),
                    AB0row = as.integer(nrow(AB0)),
                    AB0col = as.integer(ncol(AB0)),
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
    
    theta.names <- paste("theta.", subject.names, sep = "")
    alpha.beta.names <- paste(rep(c("alpha.","beta."), K),
                              rep(item.names, each = 2),
                              sep = "")
   
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)

    names <- NULL
    if(store.ability == TRUE) {
      names <- c(names, theta.names)
    }
    if (store.item == TRUE){
      names <- c(names, alpha.beta.names)
    }
    varnames(output) <- names
    attr(output,"title") <-
      "MCMCirt1d Posterior Sample"
    return(output)
    
  }
