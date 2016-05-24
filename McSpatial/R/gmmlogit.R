gmmlogit <- function(form,inst=NULL,winst=NULL,wmat=NULL,shpfile,startb=NULL,startrho=0,blockid=0,cvcrit=.0001,data=NULL,silent=FALSE) {
  library(car)

  if (identical(wmat,NULL)) {
    library(spdep)
    neighbors <- poly2nb(shpfile,queen=TRUE)
    wmat <- nb2mat(neighbors,zero.policy=TRUE)
  }
  xmat <- model.frame(form,data=data)

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
  xmat <- model.frame(form,data=data)
  y <- xmat[,1]

  xmat <- model.matrix(form,data=data)
  if (identical(inst,NULL)&identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%xmat[,-1])}
  if (identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)) {zmat <- model.matrix(inst,data=data)}
  if (!identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(model.matrix(inst,data=data), wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)) {
   cat("Warning:  list provided for inst but not winst", "\n")
   cat("inst list should include variables that are omitted from orginal explanatory variable list", "\n")
   cat("\n")
  }
  n = nrow(xmat)
  nk = ncol(xmat) 
  nk1 = nk+1

  block <- factor(blockid)
  for (i in levels(block)) {
    wmat[blockid==i,blockid!=i] <- 0
    wsum <- rowSums(wmat[blockid==i,blockid==i])
    wsum[wsum==0] <- 1
    wmat[blockid==i,blockid!=i]<- wmat[blockid==i,blockid!=i]/wsum
  }

  if (length(startb)==0) {
    logit <- glm(form,family=binomial(link="logit"),data=data)
    startb <- logit$coef
    if (silent==FALSE) {
      cat("STANDARD LOGIT ESTIMATES","\n")
      print(summary(logit))
    }
  }

  b0 <- c(startb,startrho)

  sublogit <- function(b) {
    bmat <- b[1:nk]
    rho <- b[nk1]
    xb <- xmat%*%bmat

    xstarb <- array(0,dim=n)
    u <- array(0,dim=n)
    ws <- array(0,dim=n)
    gmat <- cbind(xmat,1)
  
    for (i in levels(block)) {
      w <- wmat[block==i,block==i]
      iwmat <- solve(diag(nrow(w)) - rho*w)
      vmat <- tcrossprod(iwmat) 
      svar <- sqrt(diag(vmat))
      ws[block==i] <- diag(iwmat)/svar
      gmat[block==i,1:nk] <- iwmat%*%xmat[block==i,]/svar
      xstarb[block==i] <- gmat[block==i,1:nk]%*%bmat

      grho1 <- iwmat%*%(w%*%xstarb[block==i])
      grho2 <- diag(vmat%*%(w + t(w) - 2*rho*crossprod(w,y=NULL))%*%vmat)*xstarb[block==i]/(2*svar*svar)
      gmat[block==i,nk+1] <- grho1-as.vector(grho2)
    }
  
    p <- exp(xstarb)/(1+exp(xstarb))
    u = y-p
    du <- p*(1-p)

    for (j in seq(1:nk1)) {
      gmat[,j] <- du*gmat[,j]
      gmat[,j] <- fitted(lm(gmat[,j]~zmat))
    }

    fit <- lm(u~gmat+0)
    chmat <- fit$coef
    out <- list(chmat,gmat,ws,u,xstarb)
    names(out) <- c("chmat","gmat","ws","u","xstarb")
    return(out)
  }

  chk = cvcrit+1
  iter = 0
  chmat <- array(0,dim=nk1)
  while (chk>cvcrit|iter<=1) {
    iter = iter+1
    b0 <- b0 + chmat
    gmm <- sublogit(b0)
    chmat <- gmm$chmat
    chk = max(abs(chmat))
    print(c(iter,chk))
  }

  ws <- gmm$ws
  u <- gmm$u
  xstarb <- gmm$xstarb

  b <- b0
  gg <- solve(crossprod(gmm$gmat))
  for (j in seq(1:nk1)) {
    gmm$gmat[,j] <- abs(gmm$u)*gmm$gmat[,j]
  }
  vmat <- diag(gg%*%crossprod(gmm$gmat,y=NULL)%*%gg)
  outmat <- cbind(b, sqrt(vmat), b/sqrt(vmat), 2*(1-pnorm(abs(b)/sqrt(vmat) )) )
  v <- rownames(outmat)
  v[nrow(outmat)] = "WXB"
  rownames(outmat) <- v
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  if (silent==FALSE) {
    cat("SPATIAL GMM LOGIT ESTIMATES","\n")
    print(outmat)
  }

  out <- list(b,sqrt(vmat))
  names(out) <- c("coef","se")
  return(out)
}

  

