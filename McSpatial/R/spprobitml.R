spprobitml <- function(form,wmat=NULL,shpfile=NULL,blockid=NULL,minblock=NULL,maxblock=NULL,stdprobit=TRUE,data=NULL) {
  library(spdep)

  if (identical(data,NULL)) {
    data <- model.frame(form)
  }
  n = nrow(data)
  xmat <- model.frame(form,data=data)
  nk = ncol(xmat)+1
  if (identical(minblock,NULL)) {minblock = nk}
  if (identical(maxblock,NULL)) {maxblock = n+1}
  y <- xmat[,1]

  if (identical(blockid,NULL)) {blockid <- array(1,dim=n)}
  block <- factor(blockid)
  lblock <- levels(block)
  nblock = length(lblock)
  if (nblock>1){cat("Block diagonal W matrix will be created from shpfile; wmat will be ignored if specified","\n")}
  needw = (!identical(shpfile,NULL)&identical(wmat,NULL))|nblock>1
  if (needw==TRUE&identical(shpfile,NULL)){cat("ERROR:  shpfile required for estimation","\n")}
  tblock <- table(block)
  nbad = sum(tblock<minblock|tblock>maxblock)
  if (nbad>0) {
    badblock <- lblock[tblock<minblock|tblock>maxblock]
    lblock <- lblock[!lblock%in%badblock]
    cat("Some blocks have fewer observations than minblock","\n")
    cat("The following blocks are removed from the data set prior to estimation:","\n")
    print(tblock[rownames(tblock)%in%badblock])

    sampvar <- block%in%lblock
    data <- data[sampvar,]
    block <- block[sampvar]
    shpfile <- shpfile[sampvar,]
  }
  n = nrow(data)
  nblock = length(lblock)


  if (stdprobit==TRUE) {
    fit <- glm(form,family=binomial(link="probit"),data=data)
    cat("Standard Probit Estimates","\n")
    print(summary(fit))
  }

  makevar <- function(rho) {
    v <- array(1,dim=n)
    xstar <- model.matrix(form,data=data)

    for (i in lblock) {
      xmat <- model.matrix(form,data=data[block==i,])
      if (needw==TRUE) {
        neighbors <- poly2nb(shpfile[block==i,],queen=TRUE)
        wmat <- nb2mat(neighbors,zero.policy=TRUE)
      }
      vmat <- solve(diag(nrow(xmat)) - rho*wmat)
      xstar[block==i,] <- vmat%*%xmat
      vmat <- tcrossprod(vmat)
      v[block==i] <- sqrt(diag(vmat))
    }
    out <- list(xstar,v)
    names(out) <- c("xstar","v")
    return(out)
  }

  probitrho <- function(rho) {
    fit <- makevar(rho)
    xstar <- as.matrix(as.data.frame(fit$xstar)/fit$v)
    fit <- glm(y~xstar+0,family=binomial(link="probit"))
    xb <- fit$linear.predictors
    lvar <- sum(ifelse(y==1, log(pnorm(xb)), log(1-pnorm(xb))))
    out <- list(fit$coef,lvar)
    names(out) <- c("coef","logl")
    return(out)
  }

  logl <- function(rho) {-probitrho(rho)$logl}  
  rho = optimize(logl,lower=-.99,upper=.99)$minimum
  fit <- probitrho(rho)
  bvect <- c(fit$coef, rho)
  names(bvect) <- c(colnames(xmat),"rho")
   
  nk = length(bvect)-1
  rho = bvect[nk+1]
  b <- bvect[1:nk]
  fit <- makevar(rho)
  xb <- as.numeric(fit$xstar%*%b)/fit$v
  p <- pnorm(xb)
  g <- (dnorm(xb)^2)/(p*(1-p))
  gmat <- as.matrix((sqrt(g)/fit$v)*data.frame(fit$xstar))
  vmat1 <- solve(crossprod(gmat))

# Use numeric derivatives to calculate unconditional standard errors
  logl <- function(rho) {
    fit <- makevar(rho)
    sv <- fit$v
    xb <- as.numeric(as.matrix(fit$xstar)%*%b)/sv
    lvar <- ifelse(y==1, log(pnorm(xb)), log(pnorm(-xb)))
    return(lvar)
  }
  lvar <- sum(logl(rho))
  g <- dnorm(xb)*(y-p)/(p*(1-p))
  gmat <- as.matrix(data.frame(fit$xstar)*(g/fit$v))
  g1 <- logl(rho+.001)
  g0 <- logl(rho-.001)
  g <- (g1-g0)/.002
  vmat2 <- solve(crossprod(cbind(gmat,g)))

  semat1 <- sqrt(diag(vmat1))
  semat2 <- sqrt(diag(vmat2))
  cat("Conditional on rho","\n")
  cat("rho = ", rho, "\n")
  outmat <- cbind(bvect[1:nk], semat1, bvect[1:nk]/semat1, 2*(1-pnorm(abs(bvect[1:nk])/semat1)) )
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(outmat) <- colnames(xmat)
  print(outmat) 

  cat("Unconditional Standard Errors","\n")
  outmat <- cbind(bvect, semat2, bvect/semat2, 2*(1-pnorm(abs(bvect)/semat2)) )
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(outmat) <- c(colnames(xmat),"rho")
  print(outmat)
  cat("Number of observations = ", n, "\n")

  out <- list(bvect,lvar,vmat1,vmat2)
  names(out) <- c("coef","logl","vmat1","vmat2")
  return(out)

}



