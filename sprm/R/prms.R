prms <-
function (formula, data, a, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, center = "median", scale = "qn", usesvd = FALSE, numit=100, prec=0.01) 
  
  # PRMS Partial Robust M regression
  # (A Modified PRM R implementation compared to the "chemometrics" package).   
  # Inputs: 
  #   formula: an lm-style formula object specifying which relationship to estimate
  #   data: the data as a data frame or a list in which the y and X matrices are separate list items
  #   a: the number of PRM components to be estimated 
  #   fun (optional): the internal weighting function. Choices are "Hampel" (preferred), "Huber" or "Fair". 
  #   probp1 (optional): the 1-alpha value at which to set the first outlier cutoff
  #   probp2, probp3 (optional, only necessary if fun=="Hampel"): the 1-alpha values for cutoffs 2 and 3
  #   center (optional): how to center the data. A string that matches the R function to be used for centering 
#   scale (optional): how to scale the data. Choices are "no" (no scaling), or a string matching the R function to be used for scaling.
#   usesvd (optional, logical): whether or not to use internal singular value decomposition data compression (preferred in case of many variables). 
# Output: A list object of class "prm", having as list elements: 
#   coefficients (has the function call details as an attribute)
#   intercept
#   wy: the PRM Y space weights 
#   wt: the PRM score space weights 
#   w: the PRM combined weights 
#   scores
#   loadings 
#   fitted.values 
#   YMeans 
#   XMeans 
#   Yscales 
#   Xscales 
#   YVar: The Y pct of explained variance (overall) 
#   XVar: The X pct of explained variance (per component)
#   Written by S. Serneels, BASF SE, GVM/S, July 2013. 
#   Modified by S. Serneels, BASF Corp, GVM/T, Dec 2013.

{
#  require(vegan)
#  source("daprpr.r")
#  source("nipls.r")
  
  if(!class(formula)=="formula"){formula <- formula(formula)}
  if(is.data.frame(data) | is.list(data)){
    mt <- terms(formula, data=data)
    yname <- dimnames(attr(mt,"factors"))[[1]][1]
	if(is.list(data)){
		datnames <- names(data)
	} else {
		datnames <- colnames(data)
	}
    ic <- attr(mt, "intercept")
    if (ic==0){
      data <- tryCatch({data <- cbind(data[[which(datnames==yname)]], model.matrix(mt, data))},
                       error=function(err){
                         error <- TRUE
                         return(error)
                       }) 
    } else{
      data <- tryCatch({data <- cbind(data[[which(datnames==yname)]],model.matrix(mt, data)[,-1])},
                       error=function(err){
                         error <- TRUE
                         return(error)
                       }) 
    }
    if (is.logical(data)){
      stop("Data cannot be matched with formula.")
    } else {
      colnames(data)[1] <- dimnames(attr(mt,"factors"))[[1]][1]
    }    
  } else {
    stop("Wrong data fromat.")
  }
  
  data <- as.matrix(data)
  n <- nrow(data)
  q <- ncol(data)
  rnames <- rownames(data) # restore original rownames in the end
  rownames(data) <- 1:n  # 1:n are indices and names of w etc.
  p <- q - 1 
  
  if(length(a)>1){
    warning("Only the first element of a is used.")
    a <- a[1]
  }
  if(a>n|a>p){
    stop("The number of components a is too large.")
  }
  if (a<=0){
    stop("The number of components a has to be positive.")
  }
  if(!any(fun == c("Hampel", "Huber", "Fair"))){
    stop("Invalid weighting function. Choose Hampel, Huber or Fair for parameter fun.")
  }
  if(probp1>1|probp1<=0){
    stop("Parameter probp1 is a probability. Choose a value between 0 and 1")
  }
  if(fun=="Hampel"){
    if (!(probp1<hampelp2 & hampelp2<hampelp3 & hampelp3<=1)){
      stop("Wrong choise of parameters for Hampel function. Use 0<probp1<hampelp2<hampelp3<=1")
    }
  }
  
  if (usesvd == TRUE) {
    if (p > n) {
      dimensions <- 1
      dimension <- p - n
      cX <- colnames(data)
      ressvd <- svd(t(data[,2:q]))
      datam <- as.data.frame(cbind(data[,1],ressvd$v %*% diag(ressvd$d)))
      colnames(datam)[1] <- "Y" 
      colnames(datam)[2:(n+1)] <- paste("Xsvd",1:n,sep="")
      formula <- formula(paste("Y~",paste(colnames(datam[,2:(n+1)]),collapse="+")))
    } else {
      dimensions <- 0
      datam <- data
      warning("SVD was not used, because p doesn't exceed n")
    }
  } else {
    dimensions <- 0
    datam <- data
  }
  datamc <- daprpr(datam,center,scale)
  datac <- attr(datamc,"Center")
  datas <- attr(datamc,"Scale")
  attr(datac,"Type") <- center
  y0 <- datam[,1]
  ys <- datamc[,1]
  ns <-nrow(datamc)
  qs <- ncol(datamc)
  ps <- qs - 1
  zerows <- vector(length=0)
  wx <- sqrt(apply(datamc[,2:qs]^2, 1, sum))
  wx <- wx/median(wx)
  wy <- abs(datamc[,1])
  ###PF start modify### 
  ### wy <- wy/(1.4826*median(wy))
  if (length(wy)/2>sum(wy==0)){ # not too many zeros
    wy <- wy/(1.4826*median(wy))
  } else{
    wy <- wy/(1.4826*median(wy[wy!=0]))
  }
  ###PF end modify### 
  if(fun=="Fair"){
    wx <- 1/(1 + abs(wx/(probct*2)))
    wy <- 1/(1 + abs(wy/(probct*2)))
  }
  probct <- qnorm(probp1)
  if(fun =="Huber") {
    wx[which(wx <= probct)] <- 1
    wx[which(wx > probct)] <- probct/abs(wx[which(wx > probct)])
    wy[which(wy <= probct)] <- 1
    wy[which(wy > probct)] <- probct/abs(wy[which(wy > probct)])
  }
  if(fun =="Hampel") {
    hampelb <- qnorm(hampelp2)
    hampelr <- qnorm(hampelp3)
    wx[which(wx <= probct)] <- 1 
    wx[which(wx > probct & wx <= hampelb)] <- probct/abs(wx[which(wx > probct & wx <= hampelb)])
    wx[which(wx > hampelb & wx <= hampelr)] <- probct*(hampelr-abs(wx[which(wx > hampelb & wx <= hampelr)]))/(hampelr -hampelb)*1/abs(wx[which(wx > hampelb & wx <= hampelr)])
    wx[which(wx > hampelr)] <- 0
    wy[which(wy <= probct)] <- 1 
    wy[which(wy > probct & wy <= hampelb)] <- probct/abs(wy[which(wy > probct & wy <= hampelb)])
    wy[which(wy > hampelb & wy <= hampelr)] <- probct*(hampelr-abs(wy[which(wy > hampelb & wy <= hampelr)]))/(hampelr -hampelb)*1/abs(wy[which(wy > hampelb & wy <= hampelr)])
    wy[which(wy > hampelr)] <- 0 
  }                                                                                                                
  
  w <- wx * wy
  if(any(w<1e-6)){
    w0 <- which(w<1e-6)
    w <- replace(w,list=w0,values=1e-6)
    we <- w
  } else {
    wxe <- wx
    wye <- wy
    we <- w
  }
  dataw <- as.data.frame(datamc * sqrt(we))
  loops <- 1
  rold <- 10^-5
  difference <- 1
  while ((difference > prec) && loops < numit) {    
    res.nipls <- nipls(data=dataw,a=a)
    yp <- fitted(res.nipls)
    r <- datamc[,1] - yp
    b <- coef(res.nipls)
    Tpls <- res.nipls$scores/sqrt(we)
        if (length(r)/2>sum(r==0)){ 
      r <- abs(r)/(1.4826*median(abs(r))) 
    } else{
      r <- abs(r)/(1.4826*median(abs(r[r!=0])))
    }

    scalet = scale
    if(scale=="no"){scalet="qn"}
    dt <- daprpr(Tpls,center,scalet)
    wtn <- sqrt(apply(dt^2, 1, sum))
    wtn <- wtn/median(wtn)
    if(fun=="Fair"){
      wte <- 1/(1 + abs(wtn/(qchisq(probp1,a)*2)))
      wye <- 1/(1 + abs(r/(probct*2)))
      wye <- as.numeric(wye)
    } 
    if(fun=="Huber") {
      wte <- wtn
      wte[which(wtn <= qchisq(probp1,a))] <- 1
      wte[which(wtn > qchisq(probp1,a))] <- qchisq(probp1,a)/abs(wtn[which(wtn > qchisq(probp1,a))])
      wye <- r
      wye[which(r <= probct)] <- 1
      wye[which(r > probct)] <- probct/abs(r[which(r > probct)])
      wte[which(wtn<wt)] <- wtn[which(wtn<wt)] # ? warum ?
      wye[which(r<wye)] <- r[which(r<wye)]
      wye <- as.numeric(wye)
    }
    if(fun=="Hampel") {
      probct <- qnorm(probp1)
      hampelb <- qnorm(hampelp2)
      hampelr <- qnorm(hampelp3)
      wye <- r
      wye[which(r <= probct)] <- 1
      wye[which(r > probct & r <= hampelb)] <- probct/abs(r[which(r > probct & r <= hampelb)])
      wye[which(r > hampelb & r <= hampelr)] <- probct*(hampelr-abs(r[which(r > hampelb & r <= hampelr)]))/(hampelr -hampelb)*1/abs(r[which(r > hampelb & r <= hampelr)])
      wye[which(r > hampelr)] <- 0
      wye <- as.numeric(wye)
      
      probct <- qchisq(probp1,a)
      hampelb <- qchisq(hampelp2, a)
      hampelr <- qchisq(hampelp3, a)
      wte <- wtn
      wte[which(wtn <= probct)] <- 1 
      wte[which(wtn > probct & wtn <= hampelb)] <- probct/abs(wtn[which(wtn > probct & wtn <= hampelb)])
      wte[which(wtn > hampelb & wtn <= hampelr)] <- probct*(hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)]))/(hampelr -hampelb)*1/abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
      wte[which(wtn > hampelr)] <- 0
    }
    
    difference <- abs(sum(b^2) - rold)/rold
    rold <- sum(b^2)
    we <- wye * wte
    if(any(we<1e-6)){
      w0 <- which(we<1e-6)
      we <- replace(we,list=w0,values=1e-6)
	  zerows <- unique(c(zerows,as.numeric(names(w0))))
    }
    
    if(length(zerows)>=(n/2)){
      break
    }
    dataw <- as.data.frame(datamc * sqrt(we))
    loops <- loops + 1
  }
  if (difference > prec){
    warning(paste("Method did not converge. The scaled difference between norms of the coefficient vectors is ", round(difference, digits=4)))
  }
   
  w <- we
  w[zerows] <- 0 
  #w[setdiff(1:n,zerows)] <- we 
  #wt <- 1:n
  wt <- wte
  wt[zerows] <- 0 
  #wt[setdiff(1:length(wt),zerows)] <- wte 
  wy <- wye
  wy[zerows] <- 0 
  #wy[setdiff(1:length(wy),zerows)] <- wye 

  P <- res.nipls$loadings
  W <- res.nipls$W
  R <- res.nipls$R
  #qs <- ncol(datam)
  Tpls <- scale(datam[,2:qs],center=datac[2:qs],scale=datas[2:qs]) %*% R 
  
  if (usesvd == TRUE & dimensions == 1) {
    b <- ressvd$u %*% b 
    rownames(b) <- cX[-1]
    P <- ressvd$u %*% P 
    W <- ressvd$u %*% W
  } 

  X0 <- data[,2:q]
  Xs <- daprpr(X0,center,scale)
  Xss <- attr(Xs, "Scale")
  coef <- datas[1]/Xss*b
  if(center=="mean"){
    intercept <- mean(data[,1] - data[,2:q]%*%coef)
  } else {
    intercept <- median(data[,1] - data[,2:q]%*%coef)
  }
  
  if(!scale=="no"){
    if (center=="mean"){
      b0 <- mean(scale(data[,1],center=datac[1],scale=datas[1])-Xs%*%b)
    } else {
      b0 <- median(scale(data[,1],center=datac[1],scale=datas[1])-Xs%*%b)
    }
  } else { 
    if (center == "mean") {
      b0 <- mean(data[,1] - as.matrix(X0) %*% b)
    } else {
      b0 <- median(data[,1] - as.matrix(X0) %*% b)
    }
  }
  yfit <- as.vector(data[,2:q] %*% coef + intercept)
  resid <- as.vector(data[,1] - yfit)
  constants <- paste("cutoff1 =",probp1)
  cutoff <- probp1
  if(fun == "Hampel"){
    constants <- c(constants, paste("cutoff2 =",hampelp2), paste("cutoff3 =",hampelp3))
    cutoff <- c(cutoff, hampelp2, hampelp3)
  }
  
  dimnames(X0)[[1]] <- rnames
  dimnames(Xs)[[1]] <- rnames
  names(ys) <- rnames
  names(y0) <- rnames
  
  names(wy) <- rnames
  names(wt) <- rnames
  names(w) <- rnames
  dimnames(Tpls)[[1]] <- rnames
  dimnames(Tpls)[[2]] <- paste0("Comp", 1:(dim(Tpls)[2]))
  dimnames(W)[[2]] <- paste0("Comp", 1:(dim(W)[2]))
  dimnames(P)[[2]] <- paste0("Comp", 1:(dim(P)[2]))
  names(yfit) <- rnames
  names(resid) <- rnames
  
  inputs <- list(formula=formula, a=a,fun=fun,constants=cutoff,X0=X0,Xs=Xs,y0=y0,ys=ys,center=center,scale=scale, usesvd=usesvd)
  attr(b,"Call") <- c("PRM Regression", paste(a, "component(s)"), fun, constants, paste(center,"centering"), paste(scale,"scaling"))
  attr(coef,"Call") <- c("PRM Regression", paste(a, "component(s)"), fun, constants)
  output <- list(coefficients = coef, intercept = intercept, wy = wy, wt = wt, w = w, scores = Tpls, R=R, #weighting.vectors=W,
                 loadings = P,  fitted.values = yfit, residuals=resid, coefficients.scaled=b, intercept.scaled=b0, YMeans = datac[1], XMeans = datac[2:qs], Xscales=datas[2:qs], Yscales = datas[1], Yvar = as.vector(res.nipls$Yev), Xvar=as.vector(res.nipls$Xev),inputs=inputs)
  
  class(output) <- "prm"
  return(output)
}