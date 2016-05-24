NoncompMAR <- function(formula, Z, D, treat = NULL, data = parent.frame(),
                       n.draws = 5000, insample = TRUE, param = TRUE,
                       p.mean.c = 0, p.var.c = 100, p.mean.o = 0,
                       p.var.o = 100, smooth = 1, tie = 0.0001, mda = TRUE,
                       coef.start.c = 0, coef.start.o = 0, burnin = 0,
                       thin = 0, verbose = TRUE) {  

  ## getting the data
  call <- match.call()
  mf <- model.frame(formula, data=data, na.action='na.pass')
  X <- model.matrix(formula, data=mf)
  Y <- as.integer(model.response(mf))
  N <- length(Y)
  Z <- eval(call$Z, envir = data)
  D <- eval(call$D, envir = data)

  ## first sort by NA and then sort by treat
  treat <- eval(call$treat, envir = data)
  if (!is.null(treat)) {
    varT <- TRUE
    N11 <- length(na.omit(treat))
    indx <- 1:N
    indx <- c(indx[!is.na(treat)], indx[is.na(treat)])
    X <- X[indx,]
    Y <- Y[indx]
    Z <- Z[indx]
    D <- D[indx]
    tmp <- sort(treat, index.return = TRUE)
    difftreat <- diff(tmp$x)
    difftreat[difftreat==0] <- tie
    treat <- c(difftreat, rep(NA, N-N11))
    indx <- c(tmp$ix, (N11+1):N)
    X <- X[indx,]
    Y <- Y[indx]
    Z <- Z[indx]
    D <- D[indx]
  }
  else {
    varT <- FALSE
    N11 <- 0
  }
  Ymax <- max(na.omit(Y))
  Ymiss <- is.na(Y)
  Y[Ymiss] <- 1
  
  ## compliance status
  C <- rep(NA, N)
  C[Z==1 & D==1] <- 1
  C[Z==1 & D==0] <- 0
  C[Z==0] <- rbinom(sum(Z==0), 1, 1/2)
  res <- list(call = call, Y = Y, X = X, C = C, D = D, Z = Z,
              n.draws = n.draws)

  if (varT) {
    ## Xo = [a10 1 X Xt] where a10 for compliers without treatment
    ## The default category is the never-taker without treatment
    Xo <- cbind(0,X)
    colnames(Xo) <- c("Complier without treatment", colnames(X))
    Xo[C==1 & D==0, 1] <- 1

    ## difference matrix
    R <- rbind(c(1, rep(0, N11-2)), c(-1,1,rep(0, N11-3)))
    for (i in 1:(N11-4))
      R <- rbind(R, c(rep(0,i),-1,1,rep(0,N11-3-i))) 
    R <- rbind(R, c(rep(0, N11-3), -1, 1))
    Q <- rbind(0, diag(N11-1))
    Xt <- cbind(c(1,1,rep(0, N11-2)), Q%*%solve(R))
    Xt <- rbind(Xt, matrix(0, ncol=N11, nrow=N-N11))
    tmp <- NULL
    for (i in 1:N11)
      tmp <- c(tmp, paste("delta", i, sep=""))
    colnames(Xt) <- tmp
  }
  else {
    ## Xo = [a11 a10 1 X] where a10 for compliers without treatment
    ##                          a11 for compliers with treatment
    ## The default category is the never-taker without treatment
    Xo <- cbind(0,0,X)
    colnames(Xo) <- c("Complier with treatment", "Complier without treatment",
                      colnames(X))
    Xo[C==1 & D==1, 1] <- 1
    Xo[C==1 & D==0, 2] <- 1
  }
  ## dimension
  res$Xo <- Xo
  ncov <- ncol(X)
  ncovo <- ncovX <- ncol(Xo)
    
  ## starting values
  if(length(coef.start.c) != ncov)
    coef.start.c <- rep(coef.start.c, ncov)
  if(length(coef.start.o) != ncovo)
    coef.start.o <- rep(coef.start.o, ncovo)
 
  ## prior
  if(length(p.mean.c) != ncov)
    p.mean.c <- rep(p.mean.c, ncov)
  if(length(p.mean.o) != ncovo)
    p.mean.o <- rep(p.mean.o, ncovo)
  if(!is.matrix(p.var.c))
    p.var.c <- diag(p.var.c, ncov)
  if(!is.matrix(p.var.o))
    p.var.o <- diag(p.var.o, ncovo)
  ## prior for smooth terms
  ## putting the data and smooth stuff together
  if(varT) {
    smooth <- diag(smooth*difftreat)
    Xo <- cbind(Xo, Xt)
    ncovo <- ncol(Xo)
    p.mean.o <- c(p.mean.o, rep(0, ncovo-length(p.mean.o)))
    coef.start.o <- c(coef.start.o, rep(0, ncovo-length(coef.start.o)))
    p.var.o <- rbind(cbind(p.var.o, matrix(0, ncol = ncol(smooth),
                                           nrow = nrow(p.var.o))),
                     cbind(matrix(0, ncol = ncol(p.var.o), nrow =
                                  nrow(smooth)), smooth))
  }
  
  ## checking thinnig and burnin intervals
  if (n.draws <= 0)
    stop("`n.draws' should be a positive integer.")
  if (burnin < 0 || burnin >= n.draws)
    stop("`burnin' should be a non-negative integer less than `n.draws'.")
  if (thin < 0 || thin >= n.draws)
    stop("`thin' should be a non-negative integer less than `n.draws'.")
  keep <- thin + 1
  
  ## calling C function
  if (Ymax == 1) # binary probit
    if (param)
      allpar <- 6+ncov+ncovo
    else
      allpar <- 6
  else # ordered probit
    if (param)
      allpar <- 2*(Ymax+1)+1+ncov+ncovo+Ymax
    else
      allpar <- 2*(Ymax+1)+1
  par <- .C("MARprobit",
            as.integer(Y), as.integer(Ymiss), as.integer(Ymax),
            as.integer(Z), as.integer(D), as.integer(C), 
            as.double(X), as.double(Xo),
            as.double(coef.start.c), as.double(coef.start.o),
            as.integer(N), as.integer(n.draws),
            as.integer(ncov), as.integer(ncovo), as.integer(ncovX),
            as.integer(N11),
            as.double(p.mean.c), as.double(p.mean.o),
            as.double(solve(p.var.c)), as.double(solve(p.var.o)),
            as.integer(insample), as.integer(varT),
            as.integer(param), as.integer(mda), as.integer(burnin),
            as.integer(keep), as.integer(verbose),
            pdStore = double(allpar*(ceiling((n.draws-burnin)/keep))),
            PACKAGE="experiment")$pdStore

  ## results
  par <- matrix(par, ncol = allpar, byrow=TRUE)
  if (Ymax == 1) { # binary probit
    if (param) {
      res$coefficientsC <- par[,7:(6+ncov)]      
      if (varT){
        res$coefficientsO <- par[,(7+ncov):(6+ncov+ncovo-N11)]
        res$coefficientsS <- par[,(7+ncov+ncovo-N11):(6+ncov+ncovo)]
      }
      else
        res$coefficientsO <- par[,(7+ncov):(6+ncov+ncovo)]
    }
    res$pc <- as.matrix(par[,1])
    res$cace <- as.matrix(par[,2])
    res$itt <- as.matrix(par[,3])
    res$base <- as.matrix(par[,4:6])
    colnames(res$pc) <- "Compliance Prob."
    colnames(res$itt) <- "ITT"
    colnames(res$cace) <- "CACE"
    colnames(res$base) <- c("Baseline for compliers",
                            "Baseline for noncompliers",
                            "Baseline for all") 
  }
  else { # ordered probit 
    if (param) {
      res$coefficientsC <- par[,(2*(Ymax+1)+2):(2*(Ymax+1)+1+ncov)]
      res$coefficientsO <- par[,(2*(Ymax+1)+2+ncov):(2*(Ymax+1)+1+ncov+ncovo)]
      res$thresholds <-
        par[,(2*(Ymax+1)+2+ncov+ncovo):(2*(Ymax+1)+1+ncov+ncovo+Ymax)]
      tmp <- rep("tau", ncol(res$thresholds))
      for (i in 1:ncol(res$thresholds))
        tmp[i] <- paste("tau", i-1, sep="")
      colnames(res$thresholds) <- tmp
    }
    res$pc <- as.matrix(par[,1])
    res$cace <- as.matrix(par[,2:(Ymax+2)])
    res$itt <- as.matrix(par[,(Ymax+3):(2*(Ymax+1)+1)])
    tmp1 <- rep("ITT", ncol(res$itt))
    tmp2 <- rep("CACE", ncol(res$cace))
    for (i in 1:(Ymax+1)) {
      tmp1[i] <- paste("ITT(Y=", i-1, ")", sep="")
      tmp2[i] <- paste("CACE(Y=", i-1, ")", sep="")
    }
    colnames(res$itt) <- tmp1
    colnames(res$cace) <- tmp2
    colnames(res$pc) <- "Compliance Prob."
  }
  if (param) {
    colnames(res$coefficientsC) <- colnames(X)
    if (varT) {
      colnames(res$coefficientsO) <- colnames(Xo)[1:(ncovo-N11)]
      colnames(res$coefficientsS) <- colnames(Xt)
    }
    else
      colnames(res$coefficientsO) <- colnames(Xo)
  }
  
  class(res) <- "NoncompMAR"
  return(res)
}
