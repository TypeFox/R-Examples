## sample from the posterior distribution of a one-dimensional item
## response theory model in R using linked C++ code in Scythe.
##
## ADM and KQ 1/23/2003
## updated extensively ADM & KQ 7/28/2004
## store.ability arg added KQ 1/27/2006
## Hierarchical Subject Parameters MJM 2007-11-07
## Parameter Expansion 2008-11-18
## Grouped Subject Parameters (Wdata) started MJM,YS 2008-11

"MCMCirtHier1d" <-
  function(datamatrix, Xjdata,
           burnin = 1000, mcmc = 20000, thin=1,
           verbose = 0, seed = NA,
           theta.start = NA, a.start = NA, b.start = NA,
           beta.start=NA, b0=0, B0=.01, c0=.001, d0=.001,
           ab0=0, AB0=.25, store.item = FALSE, store.ability=TRUE,
           drop.constant.items=TRUE,
           marginal.likelihood=c("none","Chib95"), px=TRUE,
           px_a0 = 10, px_b0=10,
           ... ) {
    
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
    L <- ncol(Xjdata)       # predictors on theta Xj
    
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (J*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and call ", calling.function(), " again.\n",
              call.=FALSE)
    }
    datamatrix[is.na(datamatrix)] <- 9   
    item.names <- colnames(as.data.frame(datamatrix))
    subject.names <- rownames(as.data.frame(datamatrix))
    beta.names <- c(names(as.data.frame(Xjdata)),"sigmasq")

    ## names
    item.names <- colnames(datamatrix)
    if (is.null(item.names)){
      item.names <- paste("item", 1:K, sep="")
    }

    ## check Xj matrix and set up betastart
    Xjdata <- as.matrix(Xjdata)
    if(nrow(Xjdata) != nrow(datamatrix)) {
      cat("Error: subject covariates not of same length as datamatrix\n")
      stop("Please check data and try ",calling.function()," again.\n",call.=FALSE)
    }
    
    
    ## prior for (alpha, beta)
    holder <- form.mvn.prior(ab0, AB0, 2)
    ab0 <- holder[[1]]
    AB0 <- holder[[2]]

    ## starting values for theta error checking
    ## could use factor.score.start.check EXCEPT
    ## We have done away with eq and ineq constraints.
    if (max(is.na(theta.start))==1) { 
      theta.start <- factor.score.eigen.start(agree.mat(datamatrix), 1)
    }
    else if(is.numeric(theta.start) & length(theta.start) == J ) {
      theta.start <- theta.start * matrix(1, J, 1)  
    }
    else {
      cat("Inappropriate value of theta.start passed.\n")
      stop("Please respecify and call", calling.function(), " again.\n",
           call.=FALSE)      
    }
    ## starting values for (a, b)
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
 
    ## starting values for a and b error checking
    if (is.na(a.start)) {
      a.start <- ab.starts[,1]
    }
    else if(is.null(dim(a.start))) {
      a.start <- a.start * matrix(1,K,1)  
    }
    else if((dim(a.start)[1] != K) || (dim(a.start)[2] != 1)) {
      cat("Error: Starting value for a not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
           call.=FALSE)
    }      
    if (is.na(b.start)) {
      b.start <- ab.starts[,2]
    }
    else if(is.null(dim(b.start))) {
      b.start <- b.start * matrix(1,K,1)  
    }
    else if((dim(b.start)[1] != K) || (dim(b.start)[2] != 1)) {
      cat("Error: Starting value for b not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
              call.=FALSE)
    }    

    cat("Generating starting values (glm) for hierarchical parameters:\n")
    ## starting values are regression of theta.start on Xj
    ## or passed vector, or same value all beta
    if (max(is.na(beta.start))==1) {         # beta.start NA
      beta.start <- coef(suppressWarnings(glm.fit(Xjdata,theta.start)))
    } else if ( length(beta.start) == L ) {  # beta.start vector
      beta.start <- matrix(beta.start,L,1)
    } else if ( length(beta.start) == 1 ) {  # beta.start scalar
      beta.start <- beta.start * matrix(1,L,1)
    } else {
      cat("Error: Starting value for beta not conformable.\n")
      stop("Please respecify and call ", calling.function(), " again.\n",
           call.=FALSE)
    }
    print(beta.start)
    ## prior for beta
    holder <- form.mvn.prior(b0, B0, L)
    b0 <- holder[[1]]
    B0 <- holder[[2]]
    check.ig.prior(c0, d0)

        ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    cat("setting up posterior holder\n" )
    
    ## define holder for posterior sample
    if(store.item == FALSE & store.ability == TRUE) {
      sample <- matrix(data=0, mcmc/thin, J+L+1)
    }
    else if (store.item == TRUE & store.ability == FALSE){
      sample <- matrix(data=0, mcmc/thin, L+1 + 2*K)
    }
    else if (store.item == TRUE & store.ability == TRUE){
      sample <- matrix(data=0, mcmc/thin, L+1 + J + 2 * K)
    }
    else{
      stop("Either store.item or store.ability should be true.\n")
    }

    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    # call C++ code to draw sample
    posterior <- .C("MCMCirtHier1d",
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
                    astartdata = as.double(a.start),
                    astartrow = as.integer(length(a.start)),
                    astartcol = as.integer(1),
                    bstartdata = as.double(b.start),
                    bstartrow = as.integer(length(b.start)),
                    bstartcol = as.integer(1),
                    ab0data = as.double(ab0),
                    ab0row = as.integer(nrow(ab0)),
                    ab0col = as.integer(ncol(ab0)),
                    AB0data = as.double(AB0),
                    AB0row = as.integer(nrow(AB0)),
                    AB0col = as.integer(ncol(AB0)),
                    Xjdata = as.double(Xjdata),
                    Xjrow = as.integer(nrow(Xjdata)),
                    Xjcol = as.integer(ncol(Xjdata)),
                    betastartdata = as.double(beta.start),
                    betastartrow = as.integer(length(beta.start)),
                    betastartcol = as.integer(1),
                    b0data = as.double(b0),
                    b0row = as.integer(length(b0)),
                    b0col = as.integer(1),
                    B0data = as.double(B0),
                    B0row = as.integer(nrow(B0)),
                    B0col = as.integer(ncol(B0)),
                    c0 = as.double(c0),
                    d0 = as.double(d0),
                    storei = as.integer(store.item),
                    storea = as.integer(store.ability),
                    logmarglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    px= as.integer(px),
                    px_a0 = as.double(px_a0),
                    px_b0 = as.double(px_b0),
                    PACKAGE="MCMCpack"
                  )

    beta.names <- paste("beta.",beta.names,sep="")
    theta.names <- paste("theta.", subject.names, sep = "")
    alpha.beta.names <- paste(rep(c("a.","b."), K),
                              rep(item.names, each = 2),
                              sep = "")
   
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)
    
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
    }
    
    names <- NULL
    if(store.ability == TRUE) {
      names <- c(names, theta.names)
    }
    if (store.item == TRUE){
      names <- c(names, alpha.beta.names)
    }
    names <- c(names,beta.names)
    
    try( varnames(output) <- names)
    attr(output,"title") <-
      "MCMCirtHier1d Posterior Sample"
    attr(output,"logmarglike") <- posterior$logmarglikeholder
    return(output)
    
  }
