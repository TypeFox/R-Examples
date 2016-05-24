splogit <- function(form,inst=NULL,winst=NULL,wmat=NULL,shpfile=NULL,blockid=NULL,minblock=NULL,maxblock=NULL,data=NULL,silent=FALSE,minp=NULL) {
  library(car)
  library(spdep)

  if (!identical(winst,NULL)){wdata <- model.matrix(winst,data=data)}

  wvar = FALSE
  if (identical(wmat,NULL)&identical(shpfile,NULL)){wvar = TRUE}

  if (identical(data,NULL)) {
    data <- model.frame(form)
    xnames <- names(data)
    if (!identical(inst,NULL)){
      data1 <- model.frame(inst)
      xnames1 <- names(data1)
      newnames <- setdiff(xnames1,xnames)
      if (length(newnames)>0) {
        data <- cbind(data,data1[,newnames])
        names(data) <- c(xnames,newnames)
        xnames <- names(data)
      }
    }
    if (!identical(winst,NULL)){
      data1 <- model.frame(winst)
      xnames1 <- names(data1)
      newnames <- setdiff(xnames1,xnames)
      if (length(newnames)>0) {
        data <- cbind(data,data1[,newnames])
        names(data) <- c(xnames,newnames)
        xnames <- names(data)
      }
    }
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

  logit <- glm(form,family=binomial(link="logit"),data=data)
  if (silent==FALSE) {
    print(summary(logit))
    cat("STANDARD LOGIT ESTIMATES","\n")
  }
  xb <- logit$linear.predictors
    if (!identical(minp,NULL)) {
      xb <- ifelse(xb<qnorm(minp),qnorm(minp),xb)
      xb <- ifelse(xb>1-qnorm(minp),1-qnorm(minp),xb)
    }
  p <- exp(xb)/(1+exp(xb))
  u <- y-p
 
  gmat <- model.matrix(form,data=data)
  xnames <- colnames(gmat)
  gmat <- cbind(gmat,1)
  nk = ncol(gmat)

  for (i in lblock) {
    xmat <- model.matrix(form,data=data[block==i,])

    if (needw==TRUE) {
      neighbors <- poly2nb(shpfile[block==i,],queen=TRUE)
      wmat <- nb2mat(neighbors,zero.policy=TRUE)
    }

    if (identical(inst,NULL)&identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%xmat[,-1])}
    if (identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%(model.matrix(winst,data=data[block==i,])[,-1])) }
    if (!identical(inst,NULL)&identical(winst,NULL)) {zmat <- model.matrix(inst,data=data[block==i,])}
    if (!identical(inst,NULL)&!identical(winst,NULL)&wvar==FALSE) {
      zmat <- cbind(model.matrix(inst,data=data[block==i,]), wmat%*%(model.matrix(winst,data=data[block==i,])[,-1])) 
    }
    if (!identical(inst,NULL)&!identical(winst,NULL)&wvar==TRUE) {
      zmat <- cbind(model.matrix(inst,data=data[block==i,]), model.matrix(winst,data=data[block==i,])[,-1]) 
    }

    grad <- as.vector(p[block==i]*(1-p[block==i]))
    gmat[block==i,-nk] <- grad*xmat
    u[block==i] <- u[block==i] + gmat[block==i,-nk]%*%logit$coef
    if (wvar==FALSE){wxb <- grad*(wmat%*%xb[block==i])}
    if (wvar==TRUE){wxb <- grad*(wdata[block==i,]%*%logit$coef)}
    gmat[block==i,nk] <- wxb

    gmat[block==i,] <- zmat%*%(solve(crossprod(zmat))%*%(t(zmat)%*%gmat[block==i,]))
  }

  fit <- lm(u~gmat+0)
  v <- diag(hccm(fit))
  summat <- cbind(fit$coef,sqrt(v),fit$coef/sqrt(v),2*(1-pnorm(abs(fit$coef)/sqrt(v) )) )
  rownames(summat) <- c(xnames, "WXB") 
  colnames(summat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")

  if (silent==FALSE) {
    cat("LINEARIZED GMM LOGIT ESTIMATES","\n")
    print(round(summat,5))
    cat("Number of observations = ", n, "\n")
  }

  names(fit$coef) = c(xnames,"WXB")
  v <- sqrt(v)
  names(v) = names(fit$coef)
  out <- list(fit$coef,v,u,gmat)
  names(out) <- c("coef","se","u","gmat")
  return(out)


}

