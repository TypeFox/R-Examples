qregspiv <- function(form,wy=NULL,wmat=NULL,inst=NULL,winst=NULL,shpfile=NULL,tau=.5,rhomat=NULL,
  printsariv=FALSE,silent=FALSE,nboot=100,alpha=.05,data=NULL) {

  qdata <- model.frame(form,data=data)
  y <- qdata[,1]
  
# don't need W if wy provided and full instrument list provided, i.e., !identical(inst,NULL)&identical(inst,NULL)

  dontneedw <- !identical(wy,NULL)&!identical(inst,NULL)&identical(inst,NULL)
  
  if (identical(wmat,NULL)&dontneedw==FALSE) {
    if (identical(shpfile,NULL)) {stop("Shape file needed")}
    library(spdep)
    neighbors <- poly2nb(shpfile,queen=TRUE)
    wmat <- nb2mat(neighbors,zero.policy=TRUE)
  }
  if (identical(wy,NULL)) {wy <- as.numeric(wmat%*%as.matrix(y)) }

  qdata <- data.frame(qdata,wy)
  newform <- as.formula(form,env=qdata)

  xmat <- model.matrix(form,data=data)
  if (identical(inst,NULL)&identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%xmat[,-1])}
  if (identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(xmat, wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)) {zmat <- model.matrix(inst,data=data)}
  if (!identical(inst,NULL)&!identical(winst,NULL)) {zmat <- cbind(model.matrix(inst,data=data), wmat%*%(model.matrix(winst,data=data)[,-1])) }
  if (!identical(inst,NULL)&identical(winst,NULL)&(silent==FALSE)) {
   cat("Warning:  list provided for inst but not winst", "\n")
   cat("inst list should include variables that are omitted from orginal explanatory variable list", "\n")
   cat("\n")
  }
  zmat <- zmat[,-1]

  if (printsariv==TRUE&silent==FALSE) {
    fit <- lm(wy~zmat)
    qdata$wyhat <- fitted(fit)
    fit <- lm(update(newform,.~.+wyhat), data=qdata)
    xxmat <- summary(fit)$cov.unscaled
    yhat <- as.numeric(crossprod( t(cbind(xmat,wy)), fit$coef))
    sig2 = mean((y-yhat)^2)
    smat <- sig2*diag(xxmat)
    smat <- cbind(fit$coef,sqrt(smat),fit$coef/sqrt(smat),2*(1-pnorm(abs(fit$coef)/sqrt(smat) )) )
    rownames(smat) <- c(colnames(xmat),"WY")
    colnames(smat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    cat("\n")
    cat("IV Spatial AR Results:", "\n")
    print(smat)
    cat("sig2 =",sig2,"\n")
  }

  nk = ncol(xmat)+1
  wyboot <- wy
  zboot <- zmat
  bootdata <- qdata

  qriv <- function(dataset) {
    fit1 <- rq(wyboot~zboot,tau=tau)
    dataset$wyhat <- fitted(fit1)
    fit2 <- rq(update(newform,.~.+wyhat,env=dataset), tau=tau, data=dataset)
    return(fit2$coef)
  }

# Two-stage
  nrho = length(rhomat) 
  rhohat <- NULL
  if (nrho<=1) {
    bmat <- qriv(bootdata)
    bootmat <- array(0,dim=c(nboot,nk))
    colnames(bootmat) <- c(colnames(xmat),"wyhat")
    n = nrow(qdata)

    for (iboot in seq(1:nboot)) {
      bobs <- sample(seq(1:n),n,replace=TRUE)
      wyboot <- wy[bobs]
      zboot <- zmat[bobs,]
      bootdata <- qdata[bobs,]
      fit <- qriv(bootdata)
      xname <- names(fit) 
      bootmat[iboot,xname] <- fit
    }

    summat <- array(0,dim=c(nk,6))
    summat[,1] <- bmat
    summat[,2] <- apply(bootmat,2,sd)
    summat[,3] <- summat[,1]/summat[,2]
    summat[,4] <- 2*(1-pnorm(abs(summat[,3])) )
    qlo <- function(x) {quantile(x, alpha/2)}
    qhi <- function(x) {quantile(x, 1-alpha/2)}
    summat[,5] <- apply(bootmat,2,qlo)
    summat[,6] <- apply(bootmat,2,qhi)
    rownames(summat) <- c(colnames(xmat),"WY")
    colnames(summat) <- c("Coef.", "Bootstrap SE", "Bootstrap Z-values", "Pr(>|z|)", "Percentile-Lo", "Percentile-Hi")
    if (silent==FALSE) {
      cat("Kim and Muller Two-Stage Quantile Regression Results","\n")
      print(summat)
    }
  }

# IV
  if (nrho>1) {
    rhohat <- rhomat
    fit <- lm(wy~zmat)
    qdata$wyhat <- fitted(fit)
    newform <- as.formula(form,env=qdata)
    newform <- update(newform,newy~.+wyhat)
    for (i in seq(1:nrho)) {
      qdata$newy <- y - rhomat[i]*qdata$wy
      fit <- rq(newform,tau=tau,data=qdata)
      rhohat[i] = fit$coef[length(fit$coef)]
    }
    j = which(abs(rhohat)==min(abs(rhohat)))
    if (j==1|j==nrho) {cat("Warning:  rho is at an endpoint of rhomat","\n")}
    minrho = rhomat[j]
    if (silent==FALSE) {
      cat("Coefficients on instrumental variable for WY:","\n")
      print(cbind(rhomat,rhohat))
    }
    qdata$newy <- y - minrho*qdata$wy
    newform <- update(newform,.~.-wyhat)
    fit <- rq(newform,tau=tau,data=qdata)
    bmat <- c(fit$coef,minrho)
    smat <- crossprod(cbind(qdata$wyhat,xmat))
    e <- residuals(fit)
    h = 1.06*sd(e)*(length(y)^(-.2))
    fe <- ifelse(abs(e)<=h, .5/h, 0)
    phistar <- as.matrix(fe*qdata$wyhat)
    xstar <- as.matrix(as.data.frame(xmat)*fe)
    jmat <- solve(crossprod(cbind(phistar,xstar),cbind(qdata$wy,xmat)))
    smat <- sqrt(diag(tau*(1-tau)*(jmat%*%smat%*%t(jmat)) ))
    semat <- c(smat[2:length(smat)],smat[1])

    summat <- array(0,dim=c(nk,4))
    summat[,1] <- bmat
    summat[,2] <- semat
    summat[,3] <- bmat/semat
    summat[,4] <- 2*(1-pnorm(abs(bmat/semat)) )
    rownames(summat) <- c(colnames(xmat),"WY")
    colnames(summat) <- c("Coef.", "Std. Err.", "Z-Values", "Pr(>|z|)")
    if (silent==FALSE) {
      cat("Chernozhukov and Hansen IV Quantile Regression Results","\n")
      print(summat) 
    }
  }

    return(summat)

}

