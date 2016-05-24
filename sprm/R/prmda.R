prmda <-
  function (formula, data, a, fun="Hampel", probp1 = .95, hampelp2 = .975, hampelp3 = .999, probp4=0.01,
            yweights=TRUE, class=c("regfit", "lda"), prior=c(0.5,0.5), center = "median", scale = "qn",
            numit=100, prec=0.01) 
  {
    ## 14/10/10 IH
    ## Partial robust M regression for discriminant analysis
    ##
    ## uses daprpr, nipls, biweight, ldafitfun
    ##
    ## Inputs:
    ## formula .... an object of class formula.
    ## data ....... a data frame or list which contains the variables given in formula.
    ## a .......... the number of PRM components to be estimated in the model.
    ## fun ........ an internal weighting function for case weights. Choices are "Hampel" (preferred), "Huber" or "Fair".
    ## probp1 ..... the 1-alpha value at which to set the first outlier cutoff for the weighting function.
    ## hampelp2 ... the 1-alpha values for second cutoff. Only applies to fun="Hampel".
    ## hampelp3 ... the 1-alpha values for third cutoff. Only applies to fun="Hampel".
    ## yweights ... logical; if TRUE y weights are calculated.
    ## class ...... type of classification; choices are "regfit" or "lda".
    ## prior ...... vector of length 2 with proir porbabilities of the groups; only used if class="lda".
    ## center ..... type of centering of the data in form of a string that matches an R function, e.g. "mean" or "median".
    ## scale ...... type of scaling for the data in form of a string that matches an R function, e.g. "sd" or "qn" or alternatively "no" for no scaling.
    ## numit ...... the number of maximal iterations for the convergence of the coefficient estimates.
    ## prec ....... a value for the precision of estimation of the coefficients.
    
    #require(vegan)
    #require(pcaPP)
    #library(rrcov)
    
    if(!class(formula)=="formula"){formula <- formula(formula)}
    if(is.data.frame(data) | is.list(data)){
      mt <- terms(formula, data=data)
      yname <- dimnames(attr(mt,"factors"))[[1]][1]
      ic <- attr(mt, "intercept")
      if (ic==0){
        data <- tryCatch({data <- cbind(data[,which(colnames(data)==yname)], model.matrix(mt, data))},
                         error=function(err){
                           error <- TRUE
                           return(error)
                         }) 
      } else{
        data <- tryCatch({data <- cbind(data[,which(colnames(data)==yname)],model.matrix(mt, data)[,-1])},
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
    
    if (length(unique(data[,1]))!=2){
      stop("Wrong class labels. Only two factor levels or -1 and 1 are allowed.")
    }
    
    yorig <- data[,1]
    
    if (sum(unique(data[,1])%in%c(-1,1))<2){
      yfac <- as.factor(data[,1])
      levels(yfac) <- c("-1", "1")
      data[,1] <- as.numeric(as.character(yfac))
    }
    data <- as.matrix(data)
    n <- nrow(data)
    q <- ncol(data)
    rnames <- rownames(data)
    rownames(data) <- 1:n 
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
    if (class=="lda"){
      if (sum(prior)!=1|any(prior<=0)|length(prior)!=2){
        stop("Invalid prior probabilities. Choose two values between 0 and 1 with sum 1 for parameter prior.")
      }
    }
    if(!any(class == c("regfit", "lda"))){
      stop("Invalid classification method. Choose regfit or lda for parameter class.")
    }
    
    datam <- data
    
    scalet <- scale
    if(scale=="no"){scalet <- "qn"}
    
    ind1 <- which(datam[,1]==1)
    ind2 <- which(datam[,1]==-1)
    datamcX <- daprpr(datam[,-1],center,scale)
    datamcy <- daprpr(as.matrix(datam[,1]), center.type="mean", scale.type="sd")
    datac <- c(attr(datamcy,"Center"), attr(datamcX,"Center"))
    datas <- c(attr(datamcy,"Scale"), attr(datamcX,"Scale"))
    attr(datac,"Type") <- center
    names(datac) <- colnames(datam)
    names(datas) <- colnames(datam)
    
    datamc <- cbind(datamcy, datamcX)
    rm(datamcX)
    
    y0 <- datam[,1]
    ys <- datamc[,1]
    n1s <- length(ind1)
    n2s <- length(ind2)
    
    
    pcx1 <- prcomp(daprpr(datam[ind1,-1],center.type=center,scale.type=scalet))
    spc1 <- pcx1$sdev^2
    spc1 <- spc1/sum(spc1)
    relcomp1 <- which(spc1 - brokenstick(min(q-1,n1s)) <=0)[1]-1
    if(relcomp1==0){relcomp1 <- 1}
    r1 <- covMcd(as.matrix(pcx1$x[,1:relcomp1]))
    wx1 <- mahalanobis(as.matrix(pcx1$x[,1:relcomp1]),r1$center, r1$cov)
    # wx1 <- sqrt(apply(daprpr(as.matrix(pcx1$x[,1:relcomp1]),center,scalet)^2,1,sum))
    wx1 <- wx1/median(wx1) *qchisq(0.5,relcomp1)
    
    pcx2 <- prcomp(daprpr(datam[ind2,-1],center.type=center,scale.type=scalet))
    spc2 <- pcx2$sdev^2
    spc2 <- spc2/sum(spc2)
    relcomp2 <- which(spc2 - brokenstick(min(q-1,n2s)) <=0)[1]-1
    if(relcomp2==0){relcomp2 <- 1}
    r2 <- covMcd(as.matrix(pcx2$x[,1:relcomp2]))
    wx2 <- mahalanobis(as.matrix(pcx2$x[,1:relcomp2]),r2$center, r2$cov)
    # wx2 <- sqrt(apply(daprpr(as.matrix(pcx2$x[,1:relcomp2]),center,scalet)^2,1,sum))
    wx2 <- wx2/median(wx2) *qchisq(0.5,relcomp2)
    
    if(fun=="Fair"){
      wx1 <- 1/((1 + abs(wx1/(qchisq(probp1,relcomp1)*2)))^2) # mod: wx/(probct*2) statt wx/probct*2
      wx2 <- 1/((1 + abs(wx2/(qchisq(probp1,relcomp2)*2)))^2)
    } else if(fun =="Huber") {
      probct <- qchisq(probp1,relcomp1)
      wx1[which(wx1 <= probct)] <- 1
      wx1[which(wx1 > probct)] <- probct/abs(wx1[which(wx1 > probct)]) 
      
      probct <- qchisq(probp1,relcomp2)
      wx2[which(wx2 <= probct)] <- 1
      wx2[which(wx2 > probct)] <- probct/abs(wx2[which(wx2 > probct)])
    } else if(fun =="Hampel") {
      probct <- qchisq(probp1,relcomp1)
      hampelb <- qchisq(hampelp2,relcomp1)
      hampelr <- qchisq(hampelp3,relcomp1)
      wx1[which(wx1 <= probct)] <- 1 
      wx1[which(wx1 > probct & wx1 <= hampelb)] <- probct/abs(wx1[which(wx1 > probct & wx1 <= hampelb)])
      wx1[which(wx1 > hampelb & wx1 <= hampelr)] <- probct*(hampelr-abs(wx1[which(wx1 > hampelb & wx1 <= hampelr)]))/(hampelr -hampelb)*1/abs(wx1[which(wx1 > hampelb & wx1 <= hampelr)])
      wx1[which(wx1 > hampelr)] <- 0
      
      probct <- qchisq(probp1,relcomp2)
      hampelb <- qchisq(hampelp2,relcomp2)
      hampelr <- qchisq(hampelp3,relcomp2)
      wx2[which(wx2 <= probct)] <- 1 
      wx2[which(wx2 > probct & wx2 <= hampelb)] <- probct/abs(wx2[which(wx2 > probct & wx2 <= hampelb)])
      wx2[which(wx2 > hampelb & wx2 <= hampelr)] <- probct*(hampelr-abs(wx2[which(wx2 > hampelb & wx2 <= hampelr)]))/(hampelr -hampelb)*1/abs(wx2[which(wx2 > hampelb & wx2 <= hampelr)])
      wx2[which(wx2 > hampelr)] <- 0  
    }                                                                                                                
    
    w <- vector(length=n)
    w[ind2] <- wx2
    w[ind1] <- wx1
    if(any(w<1e-6)){
      w0 <- which(w<1e-6)
      w <- replace(w,list=w0,values=1e-6)
      we <- w
    } else {
      we <- w
    }
    dataw <- as.data.frame(datamc * we)
    colnames(dataw) <- colnames(datam)
    loops <- 1
    weold <- 10^-5
    difference <- 1
    wnorm <- vector(length=0)
    while ((difference > prec) && (loops < numit)) {   
      res.nipls <- nipls(data=dataw,a=a)
      # Tpls <- res.nipls$scores/we
      Tpls <- datamc[,-1]%*%res.nipls$R
      
      cov1 <- covMcd(as.matrix(Tpls[ind1,]))
      cov2 <- covMcd(as.matrix(Tpls[ind2,]))
      
      wlist <- int_weight(as.matrix(Tpls[ind1,]), as.matrix(Tpls[ind2,]), ind1, ind2, y0, fun, probp1, hampelp2, hampelp3, probp4, yweights=yweights,  center1=cov1$center, cov1=cov1$cov,center2=cov2$center,cov2=cov2$cov)
      we <- wlist$we
      
      dataw <- as.data.frame(datamc* we)
      colnames(dataw) <- colnames(datam)
      
      difference <- abs(sum(we^2) - weold)/weold
      weold <- sum(we^2)
      wnorm <- c(wnorm, difference)
      
      loops <- loops + 1
    }
    
    if (difference > prec){
      warning(paste("Method did not converge. The scaled difference between norms of case weights is ", round(difference, digits=4)))
    }
    
    w <- wlist$we
    wt <- wlist$wte
    if (yweights){
      wy <- wlist$wye
    } else {
      wy <- rep(1,n)
    }
    
    res.nipls <- nipls(data=dataw,a=a)
    
    b <- coef(res.nipls)
    P <- res.nipls$loadings
    W <- res.nipls$W
    R <- res.nipls$R
    
    qs <- ncol(datam)
    Tpls <- datamc[,-1] %*% R
    
    X0 <- data[,2:q]
    Xs <- daprpr(X0,center,scale)
    Xss <- attr(Xs, "Scale")
    
    coef <- datas[1]/Xss*b
    
    if (class=="lda"){
      #       # LDA im score Raum
      ind1 <- which(datam[,1]==1)
      ind2 <- which(datam[,1]==-1)
      mt1 <- apply(as.matrix(w[ind1]*Tpls[ind1,]),2,sum)/sum(w[ind1])
      mt2 <- apply(as.matrix(w[ind2]*Tpls[ind2,]),2,sum)/sum(w[ind2])
      
      Si <- matrix(0, ncol=a, nrow=a)
      for (i in 1:n){
        if (i %in% ind1){
          Si <- Si + w[i]*matrix(Tpls[i,]-mt1, ncol=1)%*% matrix(Tpls[i,]-mt1, nrow=1)
        } else if (i %in% ind2){
          Si <- Si + w[i]*matrix(Tpls[i,]-mt2, ncol=1)%*% matrix(Tpls[i,]-mt2, nrow=1)  
        }
      }
      covt <- 1/(sum(w)-2)*Si
      
      ldafit <- apply(Tpls, 1, ldafitfun, covt, mt1, mt2, prior)
      ldaclass <- (apply(ldafit, 2, which.max)-1.5)*(-2)
      ldamod <- list(cov=covt, m1=mt1, m2=mt2)
      
    } else {
      ldafit <- NULL
      ldaclass <- NULL
      ldamod <- NULL
    }
    
    if(center=="mean"){
      intercept <- mean(data[,1] - data[,2:q]%*%coef)
    } else {
      intercept <- median(data[,1] - data[,2:q]%*%coef)
    }
    
    if(!scale=="no"){
      if (center=="mean"){
        b0 <- mean(data[,1]-Xs%*%b)
      } else {
        b0 <- median(data[,1]-Xs%*%b)
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
    
    inputs <- list(a=a,formula=formula, fun=fun,constants=cutoff,X0=X0,Xs=Xs,y0=y0,ys=ys,center=center,scale=scale, prior=prior)
    attr(b,"Call") <- c("PRM Regression", paste(a, "component(s)"), fun, constants, paste(center,"centering"), paste(scale,"scaling"))
    attr(coef,"Call") <- c("PRM Regression", paste(a, "component(s)"), fun, constants)
    
    output <- list(scores = Tpls, R=R, loadings = P, 
                   w = w, wt=wt, wy=wy, Yvar = as.vector(res.nipls$Yev), Xvar=as.vector(res.nipls$Xev),
                   ldamod=ldamod, ldafit=ldafit, ldaclass=ldaclass,
                   coefficients = coef, intercept = intercept, residuales=resid, fitted.values = yfit,
                   coefficients.scaled=b, intercept.scaled=b0,
                   YMeans = datac[1], XMeans = datac[2:qs], Xscales=datas[2:qs], Yscales = datas[1], 
                   inputs=inputs)
    if(class=="lda"){
      class(output) <- "prmda"
    } else if (class=="regfit"){
      class(output) <- "prm"
    }
    return(output)
  }
