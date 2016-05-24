rda <- function(x, ...) 
{
  UseMethod("rda")
}



rda.formula <- function(formula, data = NULL, ...) 
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- rda.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("rda")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}




rda.default <- function(x, grouping=NULL, prior=NULL, 
                        gamma=NA, lambda=NA, regularization=c("gamma"=gamma, "lambda"=lambda),
                        crossval=TRUE, fold=10, train.fraction=0.5, 
                        estimate.error=TRUE, output=FALSE,
                        startsimplex=NULL, max.iter=100, trafo=TRUE,
                        simAnn=FALSE, schedule=2, T.start=0.1, halflife=50, zero.temp=0.01, alpha=2, K=100, ...)
{
  classify <- function(dataset, mu, sigma, pooled, gamma, lambda, g, p, n, prior)
  # mu            : matrix with group means as columns 
  # sigma         : (p x p x g)-Array of group covariances 
  # pooled        : pooled covariance matrix 
  # gamma, lambda : regularization parameters 
  # g, p, n       : numbers of classes, variables & observations
  {
    # compute likelihood of each datum for each group: 
    likelihood <- matrix(NA, nrow=n, ncol=g)
    dataset <- matrix(dataset,ncol=p)
    i <- 1
    singu <- FALSE
    while ((i <= g) & (!singu)){ # compute likelihoods group-wise
      reg.cov <- (1-lambda)*sigma[,,i] + lambda*pooled                   # shift towards pooled cov. 
      reg.cov <- if (is.matrix(reg.cov)) # cov is a matrix
                   (1-gamma)*reg.cov + gamma*mean(diag(reg.cov))*diag(p) # shift towards identity 
                 else # one variable  =>  cov is 1x1-matrix:
                   (1-gamma)*reg.cov + gamma*reg.cov                     # shift towards identity 
      # try to invert covariance matrix: 
      inv.trial <- try(solve(reg.cov), silent=TRUE)
      singu <- ((! is.matrix(inv.trial)) || (!is.finite(ldc<-log(1/det(inv.trial)))))
      if (!singu) likelihood[,i] <- prior[i]*.dmvnorm(dataset, mu=mu[,i], inv.sigma=inv.trial, logDetCov=ldc)
      i <- i+1
    }
    if (!singu) result <- apply(likelihood,1,function(x){order(x)[g]}) # return classifications
    else result <- rep(0,nrow(dataset)) # return definitely false classifications 
    return(result)
  }

  goalfunc <- function(paramvec=c(0.5,0.5)) 
  # depends on a (2-dim.) parameter vector, so it can be handled by minimization function. 
  # first element: gamma,  second element: lambda. 
  # returns mean misclassification rate for the bootstrap samples. 
  {
    error.rates <- numeric(fold)
    for (i in 1:fold) {
      siggi <- array(covariances[,,,i],c(p,p,g))
      mumu <- array(means[,,i], c(p,g))
      prediction <- classify(data[-train[,i],], 
                             mu=mumu, sigma=siggi, pooled=covpooled[,,i],
                             gamma=paramvec[1], lambda=paramvec[2], g=g, p=p, n=n-sum(train[,i]!=0), prior=prior)
      errors <- prediction != grouping[-train[,i]]
      group.rates <- tabulate(grouping[-train[,i]][errors],g) / test.freq[,i] #error rate by group
      group.rates[test.freq[,i]==0] <- (g-1)/g # conservative estimate for "empty" groups 
      error.rates[i] <- t(prior)%*%group.rates
    }
    return(mean(error.rates))
  }
  
  nelder.mead <- function(func, startsimplex, mini=NULL, maxi=NULL, link=function(x)(x),
                          a=1, b=0.5, g=3, r=0.5, e=.Machine$double.eps^0.5, 
                          max.iter=1e6, best.possible=-Inf, output=FALSE, 
                          simAnn=TRUE, schedule=1, T.start=NULL, halflife=200, zero.temp=0.01, alpha=2, K=500, degen=1)
  # minimizes the function `func' with several real parameters. 
  # Parameters: 
  #  func          :  the function to be minimized; must require a vector of (at least two) real-value-parameters 
  #  startsimplex  :  either a starting simplex or just the number of dimensions (for the real parameters) 
  #  mini, maxi    :  vectors of lower & upper bounds for the parameters 
  #  a             :  a > 0      `reflection coefficient' 
  #  b             :  0 < b < 1  `contraction coefficient' 
  #  g             :  g > 1      `expansion coefficient' 
  #  r             :  0 < r < 1  `reduction coefficient' 
  #  e             :  e > 0      `epsilon' (for stop criterion) 
  #  max.iter      :  (<= Inf) maximum number of iterations 
  #  best.possible :  best possible value; algorithm stops if reached. 
  #  output        :  flag for text output during calculation 
  # Parameters for Simulated Annealing: 
  #  simAnn        :  indicates whether Simulated Annealing should be used 
  #  T.start       :  starting temperature for Simulated Annealing 
  #  schedule      :  annealing schedule 1 or 2 (exponential or polynomial) 
  #  halflife      :  number of iterations until temperature is reduced to a half    (schedule I) 
  #  zero.temp     :  temperature at which it is set to zero                         (schedule I) 
  #  alpha         :  power of temperature reduction (linear, quadratic, cubic,...)  (schedule II) 
  #  K             :  number of iterations until temperature = 0                     (schedule II) 
  #  degen         :  number to substract or add to mini, maxi in degnerated simplices 
  {
    rmunif<-function(n,min=0,max=1)
    {
        dim<-length(min)
        erg<-matrix(0,nrow=n,ncol=dim)
        for (i in 1:dim) erg[,i]<-runif(n,min[i],max[i])
        if (n==1) erg<-as.vector(erg)
        return(erg)
    }
    restrict <- function(x) # forces x to be within its bounds (given by `mini' & `maxi') 
    {                       # (necessary for reflexion & expansion) 
      new.x <- x
      new.x[x<mini] <- mini[x<mini]
      new.x[x>maxi] <- maxi[x>maxi]
      return(new.x)
    }
    max.sample <- function(f) # if there is more than one worst point, one of these is sampled. 
    {
      i <- which(f==max(f))
      if (length(i)>1) i <- sample(i,1)
      return(i)
    }
    min.sample <- function(f)
    {
      i <- which(f==min(f))
      if (length(i)>1) i <- sample(i,1)
      return(i)
    }
    rlog <- function(n=1) # returns a "logarithmically distributed" random number 
    {                     # (actually equals an Exp(1)-distributed RV)            
      return(-log(runif(n)))
    }

    if (is.vector(startsimplex) && length(startsimplex)==1 && startsimplex>0 && startsimplex==trunc(startsimplex))
      startsimplex <- rbind(rep(-sqrt(((sqrt(1.5)-sqrt(0.5))^2)/2),startsimplex),diag(startsimplex))   
    p <- ncol(startsimplex)                  #  p = number of parameters 
    if (is.null(mini)) mini <- rep(-Inf,p)
    if (is.null(maxi)) maxi <- rep(Inf,p)
   
    minir<-mini
    maxir<-maxi
    for (i in 1:length(mini))
        {
        if ((mini[i]==-Inf) && (maxi[i]==Inf)) {minir[i]<- -degen; maxir[i]<- degen}
        else if ((mini[i]==-Inf) && (maxi[i]!=Inf)) minir[i]<- maxi[i]-degen
        else if ((mini[i]!=-Inf) && (maxi[i]==Inf)) maxir[i]<- mini[i]+degen
        }
    startsimplex <- t(apply(startsimplex,1,restrict))
    
    simplex<-startsimplex
    if (output) {
      if (simAnn) {
        if (schedule==1) 
            message("Performing Simulated Annealing with ", "'exponential'", " cooling schedule ",
                "(halflife=", halflife, ").")
        else{
            message("Performing Simulated Annealing with ", "'polynomial'", " cooling schedule ",
                "(alpha=", alpha, ").")
            message("Zero temperature after ", K, " iterations.")
        }
      }
      else message("Performing Nelder-Mead minimization.")
      message("Calculations for starting simplex...")
      if(.Platform$OS.type == "windows") flush.console()
    }
    f <- apply(simplex,1,function(x){func(link(x))})  #  vector of `func'-values, corresponding to the simplex 
    min.index <- min.sample(f)
    if (output) {
      message("best/worst value: ", min(f), "/", max(f))
      if(.Platform$OS.type == "windows") flush.console()
    }
    best.ever <- c(f[min.index], simplex[min.index,])
    close.enough <- ((!simAnn) && ((sd(f) <= e) || any(f <= best.possible)))
    if (simAnn) {
      if (is.null(T.start)) T.start <- min(max(f)-best.possible, max(f))
      temp <- T.start
      if (schedule==1) faktor <- 0.5^(1/halflife)
      if ((output) && (schedule==1)){ 
        message("Zero temperature after ", ceiling(log(zero.temp/T.start, base=faktor)), " iterations.")
        if(.Platform$OS.type == "windows") flush.console()
      }
    }
    i <- 1
    while ((!close.enough) && (i<=max.iter)) {  # START iterations 
      if ((simAnn) && (temp>0)) f.hot <- f + temp*rlog(p+1)    # the "annealing" function values 
      else f.hot <- f
      min.index <- min.sample(f.hot)
      max.index <- max.sample(f.hot)
      centroid <- colMeans(simplex[-max.index,])  # centroid of simplex except worst point 
     # REFLEXION: 
      x.reflex <- restrict(centroid + a*(centroid-simplex[max.index,]))
      if (output) {
        message("Reflexion")
        if(.Platform$OS.type == "windows") flush.console()
      }
      f.reflex <- func(link(x.reflex))
      if (f.reflex<best.ever[1]) best.ever <- c(f.reflex, x.reflex)
      if ((simAnn) && (temp>0)) f.reflex.hot <- f.reflex - temp*rlog(1)
      else f.reflex.hot <- f.reflex
     # If reflexion yielded improvement: EXPANSION. 
      if (f.reflex.hot <= f.hot[min.index]) {
        x.expand <- restrict(centroid + g*(x.reflex - centroid))
        if (output) { 
          message("Expansion")
          if(.Platform$OS.type == "windows") flush.console()
        }
        f.expand <- func(link(x.expand))
        if (f.expand<best.ever[1]) best.ever <- c(f.expand, x.expand)        
        if ((simAnn) && (temp>0)) f.expand.hot <- f.expand - temp*rlog(1)
        else f.expand.hot <- f.expand
       # If `expanded point' is best, it is saved, 
        if (f.expand.hot < f.hot[min.index]) {
          simplex[max.index,] <- x.expand
          f[max.index] <- f.expand
          f.hot[max.index] <- f.expand.hot
          min.index <- max.index
          max.index <- max.sample(f.hot)
        }
       # otherwise `reflected point' is saved. 
        else {
          simplex[max.index,] <- x.reflex
          f[max.index] <- f.reflex
          f.hot[max.index] <- f.reflex.hot
          min.index <- max.index
          max.index <- max.sample(f.hot)
        }
      }
     # If `reflected point' wasn't best, but better than/equal to any other point (except maximum), it is saved. 
      else if (any(f.reflex.hot <= f.hot[-max.index])) {
             simplex[max.index,] <- x.reflex
             f[max.index] <- f.reflex
             f.hot[max.index] <- f.reflex.hot
             max.index <- max.sample(f.hot)
           }
           else { # If `reflected point' was only better than worst point (but worse than others), it is saved as well. 
             if (f.reflex.hot <= f.hot[max.index]) {
               simplex[max.index,] <- x.reflex
               f[max.index] <- f.reflex
               f.hot[max.index] <- f.reflex.hot
             }
            # Since reflexion didn't improve: CONTRACTION. 
             x.contract <- centroid + b*(simplex[max.index,]-centroid)
             if (output) {
               message("Contraction")
               if(.Platform$OS.type == "windows") flush.console()
             }
             f.contract <- func(link(x.contract))
             if (f.contract<best.ever[1]) best.ever <- c(f.contract, x.contract)
             if ((simAnn) && (temp>0)) f.contract.hot <- f.contract - temp*rlog(1)
             else f.contract.hot <- f.contract
            # If `contracted' is better than (or equal to) maximum, the maximum is replaced 
             if (f.contract.hot <= f.hot[max.index]) {
               simplex[max.index,] <- x.contract
               f[max.index] <- f.contract
               f.hot[max.index] <- f.contract.hot
               min.index <- min.sample(f.hot)
               max.index <- max.sample(f.hot)
             }
             else {
              # If `contracted' was worse (greater) than maximum, the whole simplex is REDUCED: 
               for (j in (1:(p+1))[-min.index]) 
                 simplex[j,] <- simplex[min.index,]+r*(simplex[j,]-simplex[min.index,])
               if (output) {
                 message("Reduction")
                 if(.Platform$OS.type == "windows") flush.console()
               }
               f.reduce <- apply(simplex[-min.index,],1,function(x){func(link(x))})
               if (any(f.reduce<best.ever[1])) {
                 best.ever <- c(f.reduce[order(f.reduce)[1]], simplex[-min.index,][order(f.reduce)[1],])
               }
               if ((simAnn) && (temp>0)) f.reduce.hot <- f.reduce + rlog(p)  # RV is ADDED!
               else f.reduce.hot <- f.reduce
               min.index <- min.sample(f.hot)
             }
           }
           
      # Check if simplex is degeneraed?
      doubles <- apply(simplex, 2, duplicated)
      doublesdc <- colSums(doubles)
      if (any(doublesdc==p)){
       # If degenerated simplex replace by random vectors
        simplex[apply(doubles,1,any), apply(doubles,2,any)] <- 
            rmunif(sum(apply(doubles,1,any)), minir[apply(doubles,2,any)], maxir[apply(doubles,2,any)])
        f <- apply(simplex,1,function(x){func(link(x))})  
      }
      close.enough <-((sd(f) <= e) || any(f <= best.possible))
      if (output) {
        if (simAnn) 
            message(i, ".", " iteration; ",
                "temperature: ", as.character(signif(temp, 4)), "; ",
                "best/worst value: ", best.ever[1], "/", max(f))
        else 
            message(i, ".", " iteration; best/worst value: ",
                f[min.index], "/", f[max.index])
        if(.Platform$OS.type == "windows") flush.console()
      }
      i <- i+1
      if (simAnn) { # adjust temperature 
        if (schedule==1) {
          temp <- T.start * faktor^i
          if (temp < zero.temp) temp <- 0
        }
        else {
          if (i<K) temp <- T.start * (1-i/K)^alpha
          else temp <- 0
        }
        if (temp == 0) simAnn <- FALSE
      }
    }
    if (output) {
      if (close.enough) if (any(f <= best.possible)) message("Best possible value reached.")
                        else message("Converged.")
      else message("Stopped after ", i-1, " iterations.")
    }
    if ((close.enough) | any(f <= best.possible)) converged <- TRUE
    else converged <- FALSE
    simplex<-t(apply(simplex,1,link))
    result <- list(minimum=link(best.ever[-1]), value=best.ever[1], iter=i-1, converged=converged, 
                   epsilon=sd(f), startsimplex=startsimplex, finalsimplex=simplex)
    return(result)
  }
  
  crossval.sample <- function(grouping, fold=10)
  # returns more or less equally sized cross-validation-samples. 
  {
    grouping <- factor(grouping)
    g <- length(levels(grouping)) #number of groups
    if (fold > length(grouping)) fold <- length(grouping)
    cv.groups <- rep(0,length(grouping))
    groupsizes <- c(0,cumsum(summary(grouping)))
    numbers <- c(rep(1:fold, length(grouping) %/% fold), 
                 sample(fold, length(grouping) %% fold)) # group numbers to be assigned 
    for (lev in 1:g) {
      index <- which(grouping==factor(levels(grouping)[lev], levels=levels(grouping))) # indices of class "lev"
      cv.groups[index] <- sample(numbers[(groupsizes[lev]+1):groupsizes[lev+1]])
    }
    return(cv.groups)
  }

#                                                      #
#  Beginning of  _M_A_I_N_ _P_R_O_C_E_D_U_R_E_  (RDA)  #
#                                                      #
  data <- x
  rm(x)
  # if `grouping' vector not given, first data column is taken. 
  if (is.null(grouping)) {
    grouping <- data[,1]
    data <- data[,-1]
  }
  data <- as.matrix(data)
  stopifnot(dim(data)[1]==length(grouping))
  grouping <- factor(grouping)
  classes <- levels(grouping)
  if (!is.null(dimnames(data))) varnames <- dimnames(data)[[2]]
  else varnames <- NULL
  grouping <- as.integer(grouping)
  if (is.null(prior)) { 
    prior <- tabulate(grouping)
    prior <- prior / sum(prior) 
  } # frequencies as prior
  else if (all(prior == 1)) 
    prior <- rep(1 / length(classes), length(classes))   # uniform prior
  names(prior) <- classes
  dimnames(data) <- NULL
  g <- max(grouping)             # number of groups 
  p <- ncol(data)                # number of variables 
  n <- length(grouping)          # number of observations 
  if (all(is.finite(regularization)))  # no optimization 
    opti <- list(minimum=regularization, conv=FALSE, iter=0) 
  else { # optimization 
    bothpar <- (!any(is.finite(regularization)))
    n.i <- rep(round(train.fraction*n), fold) # number of observations in bootstrap (training-)samples 
    if (output) {
      message(" - RDA -")
      message(n, " observations of ", p, " variables in ", g, " classes,")
      if (crossval) message(fold, "-fold cross-validation.")
      else message(fold, " bootstrap samples of ", n.i[1], " observations each.")
      message("Class names: ", paste(classes[1:(length(classes)-1)], col=",", sep=""),
          classes[length(classes)])
      if(.Platform$OS.type == "windows") flush.console()
    }  
    # draw bootstrap/crossval samples (row indices): 
    train <- NULL
    test.freq <- NULL
      if (crossval) { # cross-validation 
        indi <- crossval.sample(grouping, fold)
        tabu <- length(indi)-tabulate(indi)
        train <- matrix(0, nrow=max(tabu), ncol=fold)
        for (i in 1:fold) {
          train[1:tabu[i],i] <- which(indi != i)
          test.freq <- cbind(test.freq, tabulate(grouping[-train[,i]], g))
        }
        n.i <- apply(train, 2, function(x) sum(x > 0))
      }
    else { # no cross-validation, but bootstrapping 
      for (i in 1:fold) {
        new <- NULL
        for (j in 1:g) new <- c(new, sample(which(grouping == j), 2))
        new <- c(new, sample((1:n)[-new], n.i - 2 * g))
        train <- cbind(train, new)
        test.freq <- cbind(test.freq, tabulate(grouping[-train[,i]], g))
        # each sample now contains at least 2 elements from each group. 
        # (samples = columns) 
      }
      remove("new")
    }
    # compute parameter estimates (mu & Sigma) for each group in each training sample: 
    means <- covariances <- covpooled <- NULL
    for (i in 1:fold) {
      means <- array(c(means, 
                       as.vector(t(as.matrix(aggregate(data[train[,i],], 
                            by = list(grouping[train[,i]]), mean)[,-1])))),
                     c(p, g, i))
      # (p x g x i)-Array ... each "slice" contains g group means as column vectors. 
      new.covar <- array(unlist(by(data[train[,i],], grouping[train[,i]], var)), c(p,p,g))
      # (p x p x g)-Array, each slice is covariance matrix of one group. 
      covariances <- array(c(covariances,new.covar), c(p,p,g,i))
      # (p x p x g x i)-Array, each "hyperslice" contains g covariance matrices, as above. 
      new.cp <- array(new.covar, c(p*p,g))
      weights <- (tabulate(grouping[train[,i]])-1) / (n.i[i]-g) 
      # weights proportional to fraction of group in sample 
      covpooled <- array(c(covpooled, matrix(new.cp %*% weights, p, p)), c(p, p, i))
      # (p x p x i)-Array, each slice contains pooled covariance for i-th bootstrap sample. 
    }
    remove(list=c("new.covar","new.cp","weights"))
    if (bothpar) { # optimization over both parameters 
      if (is.null(startsimplex)) {
        #startsimplex <- matrix(rbeta(6,1,1),ncol=2) # Beta-RVs for Startsimplex
        startsimplex <- cbind(c(runif(1,1/11,4/11),runif(1,4/11,7/11),runif(1,7/11,10/11)),
                              c(runif(1,1/11,4/11),runif(1,4/11,7/11),runif(1,7/11,10/11)))
        perm <- cbind(c(1,3,2), c(2,1,3), c(3,1,2), c(2,3,1))
        startsimplex[,2] <- startsimplex[perm[,sample(4,1)],2]
      }
      if (trafo) {  # use transformation in Nelder-Mead. 
        linkfunc <- function(x)
            return(1/(1+exp(-x))) # sigmoidal transformation function (for both parameters) 
        linkinverse <- function(x)
            return(-log(1/x-1))   # inverse of link function 
        #startsimplex <- matrix(rnorm(6,0,1),ncol=2) # Normal-RVs for Startsimplex 
        startsimplex <- linkinverse(startsimplex) # transform startsimplex
        mini <- c(-Inf, -Inf)
        maxi <- c(Inf, Inf)
      }
      else {        # do not use transformation. 
        linkfunc <- function(x)
            return(x) # identity function 
        #linkinverse <- linkfunc
        mini=c(0,0)
        maxi=c(1,1)
      }
      dimnames(startsimplex) <- list(NULL, c("gamma","lambda"))
      # NELDER-MEAD #
      opti <- nelder.mead(goalfunc, startsimplex, mini=mini, maxi=maxi, 
        link=linkfunc, max.iter=max.iter, best.possible=0, out=output, 
        simAnn=simAnn, schedule=schedule, T.start=T.start, 
        halflife=halflife, zero.temp=zero.temp, alpha=alpha, K=K)
    }
    else { # optimization over single parameter 
      logit <- function(x)  return(1/(1+exp(-x)))
      tryval <- logit(seq(-4,4,le=12)) # values to try first 
      if (is.na(regularization[1])) {
        if (output) {
          message("Optimizing gamma...")
          if(.Platform$OS.type == "windows") flush.console()
        }
        goalfu2 <- function(x)
            return(goalfunc(c(x, regularization[2])))
        err <- apply(matrix(tryval, ncol=1), 1, goalfu2)
        fromto <- c(0, tryval, 1)[which.min(err) + c(0,2)]
        minimize <- optimize(goalfu2, fromto)
        opti <- list(minimum=c(minimize$minimum, regularization[2]), 
                     value=minimize$objective, conv=TRUE, iter=-1)
      }
      else {
        if (output) {
          message("Optimizing lambda...")
          if(.Platform$OS.type == "windows") flush.console()
        }
        goalfu2 <- function(x)
        {return(goalfunc(c(regularization[1], x)))}
        err <- apply(matrix(tryval,ncol=1),1,goalfu2)
        fromto <- c(0,tryval,1)[which.min(err)+c(0,2)]
        minimize <- optimize(goalfu2,fromto)
        opti <- list(minimum=c(regularization[1], minimize$minimum), 
                     value=minimize$objective, conv=TRUE, iter=-2)     
      }
    }
  }
  opt.par <- opti$minimum; names(opt.par) <- c("gamma","lambda")
  if (output) {
    message("Regularization parameters:\n gamma: ", round(opt.par[1],5), 
        "  lambda:", round(opt.par[2],5))
    if(.Platform$OS.type == "windows") flush.console()
  }
  # compute parameters for complete data: 
  means <- t(as.matrix(aggregate(data,by=list(grouping),mean)[,-1]))
  dimnames(means) <- list(varnames,classes)
 # covariances <- array(unlist(by(data[train[,i],],grouping[train[,i]],var)),
 #                      c(p,p,g), dimnames=list(varnames,varnames,classes))
  covariances <- array(unlist(by(data,grouping,var)),
                       c(p,p,g), dimnames=list(varnames,varnames,classes))
 # weights <- (tabulate(grouping)-1) / (n.i-g)
  weights <- (tabulate(grouping)-1) / (n-g)
  covpooled <- matrix(array(covariances, c(p*p,g)) %*% weights, p, p,
                      dimnames = list(varnames, varnames))
  # predict training data, compute apparent error rate: 
  if (estimate.error) {
    errors <- classify(data, means, covariances, covpooled, opt.par[1], opt.par[2], g=g, p=p, n=n, prior=prior) != grouping
    group.rates <- tabulate(grouping[errors],g) / tabulate(grouping)
    APER <- as.vector(t(prior) %*% group.rates)
    if (output) 
        message("Apparent error rate (APER) for training data: ", 
            round(APER * 100, 3), "%")  
  }
  else APER <- NA
  if (crossval) err <- c("APER"=APER, "crossval"=opti$value)
  else err <- c("APER"=APER, "bootstrap"=opti$value)
  result <- list(call=match.call(), 
                 regularization=opt.par, classes=classes, prior=prior, error.rate=err,
                 varnames=varnames,
                 means=means, covariances=covariances, covpooled=covpooled, 
                 converged=opti$conv, iter=opti$iter)
  class(result) <- "rda"
  return(result)
}

predict.rda <- function(object, newdata, posterior=TRUE, aslist=TRUE, ...)
{
  classify <- function(dataset, mu, sigma, pooled, gamma, lambda, g, p, n, prior)
  # difference to `classify'-function above is that LIKELIHOODS are returned instead of classifications
  # mu            : matrix with group means as columns 
  # sigma         : (p x p x g)-Array of group covariances 
  # pooled        : pooled covariance matrix 
  # gamma, lambda : regularization parameters
  # g, p, n       : numbers of classes, variables & observations
  {
    # compute likelihood of each datum for each group: 
    likelihood <- matrix(nrow=n, ncol=g)
    for (i in 1:g){ 
      reg.cov <- (1-lambda) * sigma[,,i] + lambda * pooled               # shift towards pooled cov. 
      reg.cov <- if (is.matrix(reg.cov)) # cov is a matrix
                   (1-gamma)*reg.cov + gamma*mean(diag(reg.cov))*diag(p) # shift towards identity 
                 else # one variable  =>  cov is 1x1-matrix:
                   (1-gamma)*reg.cov + gamma*reg.cov                     # shift towards identity 
      likelihood[,i] <- prior[i] * .dmvnorm(dataset, mu = mu[,i], inv.sigma = solve(reg.cov))
    }
    return(likelihood)
  }

  p <- dim(object$means)[1]
  g <- dim(object$means)[2]
  
      if (!inherits(object, "rda")) 
        stop("object not of class", " 'rda'")
        if (!is.null(Terms <- object$terms)) {
        if (missing(newdata)) 
            newdata <- model.frame(object)
            else {
                newdata <- model.frame(as.formula(delete.response(Terms)), 
                newdata, na.action = function(x) x, xlev = object$xlevels)
            }
            x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
            xint <- match("(Intercept)", colnames(x), nomatch = 0)
            if (xint > 0) 
                x <- x[, -xint, drop = FALSE]
            }
        else {
        if (missing(newdata)) {
             if (!is.null(sub <- object$call$subset)) 
                    newdata <- eval.parent(parse(text = paste(deparse(object$call$x, 
                      backtick = TRUE), "[", deparse(sub, backtick = TRUE), 
                    ",]")))
                else newdata <- eval.parent(object$call$x)
                if (!is.null(nas <- object$call$na.action)) 
                 newdata <- eval(call(nas, newdata))
            }
            if (is.null(dim(newdata))) 
             dim(newdata) <- c(1, length(newdata))
            x <- as.matrix(newdata)
        }

  
  
  #if (!any(is.null(colnames(newdata)),is.null(object$varnames))) { # (both colnames & varnames are given) 
  #  if(all(is.element(object$varnames,colnames(newdata)))){        # (varnames is a subset of colnames)   
  #    newdata <- as.matrix(newdata[,object$varnames])
  #  }
  #}
  #if(is.vector(newdata)) newdata <- matrix(newdata, ncol=1)
  n <- dim(newdata)[1]
  newdata<-x
  likeli <- classify(newdata, object$means, object$covariances, object$covpooled,
                     object$regul[1], object$regul[2], g=g, p=p, n=dim(newdata)[1], 
                     prior=object$prior)
  colnames(likeli) <- object$classes
  classi <- apply(likeli, 1, function(x) order(x)[g])
  classi <- factor(object$classes[classi], levels = object$classes)
  postmat <- if (posterior) likeli / rowSums(likeli)
             else NULL
  if (aslist) result <- list("class" = classi, "posterior" = postmat)
  else { 
    result <- classi
    attr(result, "posterior") <- postmat
  }
  return(result)
}


print.rda <- function(x, ...)
{
  #cat(" - RDA - \n")
  cat("Call:", "\n")
  print(x$call, ...)
  cat("\nRegularization parameters:", "\n")
  #cat("gamma:", round(x$regu[1],5), " lambda:", round(x$regu[2],5), "\n")
  print(x$regu, ...)
  #cat("\nClass prior:", "\n")
  cat("\nPrior probabilities of groups:", "\n")
  print(x$prior, ...)
  cat("\nMisclassification rate:", "\n")
  cat("       apparent:",
    ifelse(is.na(x$error.rate[1]), "--", as.character(round(x$error.rate[1] * 100, 3))), "%\n")
  if (length(x$error.rate) > 1)
    cat(ifelse(names(x$error.rate)[2] == "crossval", 
        "cross-validated:", "   bootstrapped:"), 
        as.character(round(x$error.rate[2] * 100, 3)), "%\n")
  invisible(x)
}


plot.rda <- function(x, textplot=FALSE, ...)
{
  parpty <- par("pty")
  par(pty="s")
  if(textplot) {
    plot(c(0,1),c(0,1), type="n", axes=FALSE, xlab="groups", ylab="covariances", ...)
    textcol <- "darkgrey"
    textsize <- 1.5
    textshift <- 0.02
    text(0+textshift, 0+textshift, "QDA", adj=c(0,0), cex=textsize, col=textcol)
    text(1-textshift, 0+textshift, "LDA", adj=c(1,0), cex=textsize, col=textcol)
    #text(0+textshift, 1-textshift, "cond. indep.", adj=c(0,1), cex=textsize, col=textcol)
    #text(1-textshift, 1-textshift, "nearest mean", adj=c(1,1), cex=textsize, col=textcol)
    text(0.5, 1-textshift, "i.i.d. variables", adj=c(0.5,1), cex=textsize, col=textcol)
    axis(1, at=c(0,1), labels=c("unequal", "equal"))
    axis(2, at=c(0,1), labels=c("correlated", "diagonal"))
  }
  else{
    plot(c(0,1),c(0,1), type="n", axes=FALSE, xlab=expression(lambda), ylab=expression(gamma), ...)
    axis(1)
    axis(2)
  }
  lines(c(0,1,1,0,0), c(0,0,1,1,0), col="grey")
  lines(c(0,1), rep(x$regu[1],2), col="red1", lty="dotted")
  lines(rep(x$regu[2],2), c(0,1), col="red1", lty="dotted")
  points(x$regu[2], x$regu[1], pch=18, col="red2")
  par(pty=parpty)
  invisible(x$regu)
}


.dmvnorm <- function(x, mu=NA, inv.sigma=NA, logDetCov=NA)    
# Density of a Multivariate Normal Distribution 
# works for matrices as well as for vectors     
# (matrices are evaluated row-wise)             
#   !!   supply inverse of covariance   !!      
# `logDetCov' = log(det(Cov)) = log(1/det(inv.Cov)) = -1*log(det(inv.Cov))
{
  if (is.vector(x)) x <- t(x) # x is treated as 1 obsevation of (length(x)) variables
  if (is.na(logDetCov)) logDetCov <- -log(det(inv.sigma))
  singledens <- function(x, M=mu, IS=inv.sigma, ldc=logDetCov)  # density function for a single vector 
  {
    xm <- x - M
    return(as.numeric(exp(- 0.5 * 
           (length(M) * log(2*pi) 
            + ldc
            + (t(xm) %*% IS %*% xm)))))
  } 
  return(apply(x, 1, singledens))#, M=mu, IS=inv.sigma, ldc=logDetCov))
}  
