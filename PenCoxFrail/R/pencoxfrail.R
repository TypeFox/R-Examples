## roxygen2 package
##############################################################
##############################################################
##############################################################
##############################################################
pencoxfrail<- function(fix=formula, rnd=formula, vary.coef=formula, 
                       data, xi, adaptive.weights = NULL, control = list()){
  
  est <- est.pencoxfrail(fix=fix,rnd=rnd,vary.coef=vary.coef,data=data,
                         xi=xi,adaptive.weights=adaptive.weights,control=control)
  
  est$StdDev <- est$Q
  est$call <- match.call()
  class(est) <- "pencoxfrail"
  est
}



est.pencoxfrail <- function(fix, rnd, vary.coef, data, xi, adaptive.weights = NULL, control = list() )
{
  if(grepl("\\*", fix[3]))
    stop("Usage of '*' not allowed in formula! Please specify the corresponding variables separately.")  
  
  ic.dummy<-attr(terms(fix),"intercept")  
  
  Resp <- model.response(model.frame(fix, data))
  
  if(ncol(Resp)==2)
  {
    y0 <- rep(0,nrow(Resp))  
    y <- Resp[,1]
    delta <- Resp[,2]
  }else{
    y0 <- Resp[,1]  
    y <- Resp[,2]
    delta <- Resp[,3]
  }
  
  xl <- 0
  N <- length(y)
  
  vary.names <- attr(terms(vary.coef),"term.labels")
  
  ind1 <- sapply(data[,vary.names,drop=FALSE],is.factor)
  ind2 <- sapply(vary.names,factor.test)
  ind <- (ind1+ind2)>0
  
  number.cat <- sum(ind)
  number.noncat <- sum(!ind)
  
  cat.present <- number.cat>0
  noncat.present <- number.noncat>0 

  if(cat.present)
  {
    vary.coef.cat <- formula(paste("~ 1+",paste(vary.names[ind],collapse="+")))
  }else{
    vary.coef.cat <- formula("~ 1")
  }
  
  if(noncat.present)
  {
    vary.coef <- formula(paste("~ 1+",paste(vary.names[!ind],collapse="+")))
  }else{
    vary.coef <- formula("~ 1")
  }
  
  if(cat.present)
  {
    levels.vec <- numeric()
    for(i in which(ind==TRUE))
      levels.vec[i] <- nlevels(data[,vary.names[i]])-1
    cum.levels <- c(0,cumsum(levels.vec))
  }
  
  control<-do.call(pencoxfrailControl, control)
  
  if(is.null(control$xr))
  {
    xr <- max(y)
  }else{
    xr <- control$xr
  }
  
  if(ic.dummy!=1){
    fix.help <- update(fix,~ .+1)
    X <- model.matrix(fix.help, data)[,-1]
  }else{
    X <- model.matrix(fix, data)  
  }
  
  ### remove intercept
  X <- X[,-1]
  
  ## set up random formula
  rnd.len<-length(rnd)
  
  rndformula <- as.character(rnd)
  
  trmsrnd <- terms(rnd[[1]])
  
  if(!is.factor(data[,names(rnd)[1]]))
  {
    data[,names(rnd)[1]] <- as.factor(data[,names(rnd)[1]])
    warning("Cluster variable should be specified as a factor variable!")  
  }
  
  newrndfrml <- "~ -1"
  newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")
  
  if(length(attr(trmsrnd, "variables"))>1)
  {
    newrndfrml <- paste(newrndfrml,  
                        paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                          paste(lbl, names(rnd)[1], sep=":")
                        }), collapse=" + "), sep="+") 
  }
  
  W_start <- model.matrix(formula(newrndfrml), data)
  
  rnlabels <- terms(formula(newrndfrml))
  random.labels <- attr(rnlabels,"term.labels")
  ran.tab <- table(data[,colnames(data)==(names(rnd)[1])])   
  n <- length(ran.tab)
  s <- dim(W_start)[2]/n
  
  if(s>1)
  {
    W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
    for (i in 2:n)
      W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
  }else{
    W<-W_start
  }
  subject.names<-names(rnd)  

    XW <- cbind(X,W)
  
  ### time varying coef
  if(attr(terms(vary.coef), "intercept")==0)
  {
    variables<-attr(terms(vary.coef),"term.labels")
    vary.coef<- paste(rownames((attr(terms(vary.coef),"factors")))[1],"~ +1",sep="")
    for (ir in 1:length(variables))
      vary.coef <- paste(vary.coef, variables[ir], sep="+")
    vary.coef <-formula(vary.coef)
  }
  
  smooth <- control$smooth
  start <- control$start
  exact <- control$exact
  
  ridge.pen <- control$ridge.pen
  
  nbasis <- smooth$nbasis
  diff.ord <- 2
  spline.degree <- 3
  penal<-smooth$penal
  
  U <- model.matrix(vary.coef, data)
  m <- ncol(U)-1
  
  U2 <- model.matrix(vary.coef.cat, data)
  U2 <- U2[,-1]
  m2 <- ncol(U2)

  U.both <- cbind(U,U2)
  
  ## inititalize adaptive weights
  if(is.null(adaptive.weights))
    adaptive.weights <- matrix(1,nrow=((m+m2+1)*nbasis),ncol=2)
  
  ## standardization
  mean.vec <- colMeans(U.both[,2:ncol(U.both),drop=FALSE])
  sd.vec <- apply(U.both[,2:ncol(U.both),drop=FALSE],2,sd)
  
  if(control$center)
    U.both[,2:ncol(U.both)] <- scale(U.both[,2:ncol(U.both)],scale=F)
  
  if(control$standardize)
    U.both[,2:ncol(U.both)] <- scale(U.both[,2:ncol(U.both)],center=F)
  
  U.each <- t(apply(U.both,1,rep,each=nbasis))
  
  U.double.each <- matrix(NA,N,nbasis^2*ncol(U.both)^2)
  for(j in 1:(m+m2+1))
    U.double.each[,((j-1)*(m+m2+1)*nbasis^2+1):(j*(m+m2+1)*nbasis^2)] <- t(apply(t(apply(U.both*U.both[,j],1,rep,each=nbasis)),1,rep,nbasis))
  
  if(ncol(U.both)==1)
  {
    if(colnames(U.both)=="(Intercept)")
      stop("No terms to select! Use coxph or coxme!")
  }
  
  ###
  
  if(is.null(control$q_start))
  {
    control$q_start<-rep(0.1,sum(s))
    if(sum(s)>1)
      control$q_start<-diag(control$q_start)
  }
  
  Q.start<-control$q_start
  
  lin<-ncol(X)
  
  zeta <- control$zeta
  c.app <- control$c.app
  ## start with estimation
  
  if(control$print.iter)
    message("Iteration ",1)
  
  ## generate penalty matrix
  D<-diag(nbasis)
  D.1 <- diff(D)
  D<-diff(D.1);t.D<-t(D)
  
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(nbasis-3)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  knots<-knots[1:nbasis]
  
  Pen <- t.D%*%solve(D%*%t.D)
  B.unpen.fact<-cbind(1,knots)
  
  K.part <- cbind(B.unpen.fact,Pen)
  
  K.smooth <- crossprod(K.part)
  K.const <- crossprod(D.1%*%K.part)
  
  if(cat.present)
  {
    K.eucl <- list()
    for(ee in 1:number.cat)
    {
      K.list <- bdiag(replicate(levels.vec[ee],K.part,simplify = F))
      K.eucl[[ee]] <- as.matrix(crossprod(K.list))
    }
  }
  
  K.base <- c(0,0,rep(1,nbasis-2))
  K.ridge <- rep(c(0,0,rep(1,nbasis-2)),m+m2)
  
  if(!(diff.ord<spline.degree))
    stop("Order of differences must be lower than degree of B-spline polynomials!")
  
  Phi.base<-bs.design(y,xl=xl,xr=xr,spline.degree=spline.degree,nbasis=nbasis)  
  colnames(Phi.base)<-paste("baseline",rep(1:dim(Phi.base)[2],each=1), sep=".")
  
  Phi.base  <-cbind(Phi.base%*%B.unpen.fact,Phi.base%*%Pen)
  
  Phi<-matrix(0,N,(m+m2+1)*nbasis)
  for(j in 1:(m+m2+1))
    Phi[,((j-1)*nbasis+1):(j*nbasis)] <- Phi.base * U.both[,j]
  
  #######
  time.grid <- seq(0,max(y),by=exact)
  time.length <- length(time.grid)
  
  T.mat <- matrix(NA,time.length,N)
  event.ma <- cbind(y0,y)
  T.mat[1:length(time.grid),1:N] <- apply(event.ma,1,function(z,time.grid){(time.grid<z[2] & time.grid>=z[1])},time.grid=time.grid)
  
  Phi.big.temp <- bs.design(time.grid,xl=xl,xr=xr,spline.degree=spline.degree,nbasis=nbasis)
  colnames(Phi.big.temp) <- paste("baseline",rep(1:dim(Phi.base)[2],each=1), sep=".")
  Phi.big.temp  <-cbind(Phi.big.temp%*%B.unpen.fact,Phi.big.temp%*%Pen)
  
  Phi.big <- t(apply(Phi.big.temp,1,rep,m+m2+1))
  
  Phi.double.big <- matrix(NA,length(time.grid),(m+m2+1)*nbasis^2)
  for(j in 1:nbasis)
    Phi.double.big[,((j-1)*nbasis*(m+m2+1)+1):(j*nbasis*(m+m2+1))] <- t(apply(Phi.big.temp*Phi.big.temp[,j],1,rep,m+m2+1))
  Phi.double.big <- (t(apply(Phi.double.big,1,rep,m+m2+1)))
  
  ## define transposed matrices
  t.U.each <- t(U.each)
  t.XW <- t(XW)

  ## define products
  delta.Phi <- delta %*% Phi
  
  dimb <- ncol(Phi)
  dimr <- n*s
  
  if(is.null(start))
    start <- rep(0,dimb + lin + dimr)
  
  Delta <- matrix(NA,control$max.iter+1, dimb + lin + dimr)
  Delta[1,] <- start
  
  Q<-list()
  Q[[1]]<-Q.start
  
  alpha <- start[1:dimb]
  if(lin>0)
  {
    fixef <- start[(1+dimb):(dimb+lin)]
  }else{
    fixef <- NULL
  }
  ranef <- start[(1+dimb+lin):(dimb+lin+dimr)]
  alphama <- matrix(alpha,ncol(U.each),nrow(U.each))
  multi.obj <- IntegrMultiCpp(alphama, c(fixef,ranef), Phi.big,
                              t.U.each, U.each, T.mat, Phi.double.big, U.double.each, XW)
  multi.obj$int.ma <- exact*multi.obj$int.ma
  multi.obj$int.array <- exact*multi.obj$int.array
  multi.obj$eta <- drop(multi.obj$eta)
  
  int.vec.fix <-  exact * rowSums(multi.obj$help.calc)
  
  Fisher.upp <- matrix(multi.obj$int.array, dimb, dimb)
  
  P <- matrix(0,dimb+lin+dimr,dimb+lin+dimr)
  diag(P)[1:nbasis] <- penal * K.base
  
  if(noncat.present)
  {
    alpha.ma.small <- matrix(alpha[(nbasis+1):((m+1)*nbasis)],ncol=m,nrow=nbasis)
    pen.1.vec <- zeta * sqrt(nbasis-1) * apply(alpha.ma.small,2,penal.fct.inv,c.app = c.app, K = K.const) * adaptive.weights[1:m,1]
    pen.2.vec <- (1-zeta) * sqrt(nbasis) * apply(alpha.ma.small,2,penal.fct.inv,c.app = c.app, K = K.smooth) * adaptive.weights[1:m,2]
    for(r in 1:m)
    {
      P[(r*nbasis+1):((r+1)*nbasis),(r*nbasis+1):((r+1)*nbasis)] <-
        xi * ( pen.1.vec[r] * K.const + pen.2.vec[r] * K.smooth) 
    }
  }
  
  if(cat.present)
  {
    levels.ind <- rep(1:number.cat,levels.vec)
    alpha.ma2 <- matrix(alpha[((m+1)*nbasis+1):((m+m2+1)*nbasis)],nbasis,m2)
    eucl.norm.vec <- rep(0,number.cat)
    for(k in 1:number.cat)
      eucl.norm.vec[k] <- sqrt(levels.vec[k]*nbasis)*penal.fct.inv(as.vector(alpha.ma2[,(levels.ind==k)]),K.eucl[[k]],c.app)
    
    pen.1.vec <- zeta * sqrt(nbasis-1) * apply(alpha.ma2,2,penal.fct.inv,c.app = c.app, K = K.const) * adaptive.weights[(m+1):(m+m2),1]
    pen.2.vec <- (1-zeta) * eucl.norm.vec * adaptive.weights[(m+cumsum(levels.vec)),2] 
    
    for(r in 1:m2)
    {
      P[((r+m)*nbasis+1):((r+m+1)*nbasis),((r+m)*nbasis+1):((r+m+1)*nbasis)] <-
        xi * ( pen.1.vec[r] * K.const)
    }
    
    for(r in 1:number.cat)
      P[((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis),((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis)]  <- P[((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis),((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis)] + xi * ( pen.2.vec[r] * K.eucl[[r]])
  }
  
  if(s==1)
  {
    diag(P)[(dimb+lin+1):(dimb+lin+dimr)] <- 1/Q.start
  }else{
    Q_inv.start <- chol2inv(chol(Q.start))
    P <- matrix(0,dimb+lin+dimr,dimb+lin+dimr)
    for(jf in 1:n)
      P[(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s),(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s)] <- Q_inv.start
  }
  
  diag(P)[(nbasis+1):((m+m2+1)*nbasis)] <- diag(P)[(nbasis+1):((m+m2+1)*nbasis)] + ridge.pen * K.ridge
  
  score <- c(delta.Phi - multi.obj$t.eta %*% multi.obj$int.ma,
             t.XW %*%(delta -  (int.vec.fix * multi.obj$eta))) - P %*% start
  
  Fisher <- matrix(NA,dimb+lin+dimr,dimb+lin+dimr)
  
  Fisher[1:dimb,1:dimb] <- -Fisher.upp
  
  Fisher[(dimb+1):(dimb+lin+dimr),1:dimb] <-  - t.XW %*% (multi.obj$int.ma*multi.obj$eta)
  
  Fisher[1:dimb,(dimb+1):(dimb+lin+dimr)] <- 
    t(Fisher[(dimb+1):(dimb+lin+dimr),1:dimb])
  
  Fisher[(dimb+1):(dimb+lin+dimr),(dimb+1):(dimb+lin+dimr)] <- 
    -t.XW %*% (XW*multi.obj$eta*int.vec.fix)
  
  ###
  
  InvFisher<-try(chol2inv(chol(P - Fisher)),silent=T)
  if(class(InvFisher)=="try-error")
    InvFisher<-try(solve(P - Fisher),silent=T)
  if(class(InvFisher)=="try-error" || sum(is.na(InvFisher))>0)
    stop("Current Fisher matrix not invertible")  
  
  Delta[2, ] <- start + InvFisher %*% score
  
  ### random effects variance
  alpha <- Delta[2,1:dimb]
  if(lin>0)
  {
    fixef <- Delta[2,(dimb+1):(dimb+lin)]
  }else{
    fixef <- NULL
  }
  ranef <- Delta[2,(1+dimb+lin):(dimb+lin+dimr)]
  
  if(s==1)
  {
    Q1 <- sum(diag(InvFisher)[(dimb+1):(dimb+lin+dimr)])+sum(ranef^2)
  }else{ 
    Q1<-InvFisher[(dimb+lin+1):(dimb+lin+s),(dimb+lin+1):(dimb+lin+s)]+ranef[1:s]%*%t(ranef[1:s])
    for (i in 2:n)
      Q1<-Q1+InvFisher[(dimb+lin+(i-1)*s+1):(dimb+lin+i*s),(dimb+lin+(i-1)*s+1):(dimb+lin+i*s)]+ranef[((i-1)*s+1):(i*s)]%*%t(ranef[((i-1)*s+1):(i*s)])
  }  
  Q1<-Q1/n
  
  Q[[2]]<-Q1
  ########################################################################
  ########################################################################
  ########################################################################
  for(ll in 2:control$max.iter)
  {
    if(control$print.iter)
      message("Iteration ",ll)
    
    alphama <- matrix(alpha,ncol(U.each),nrow(U.each))
    multi.obj <- IntegrMultiCpp(alphama, c(fixef,ranef), Phi.big,
                                t.U.each, U.each, T.mat, Phi.double.big, U.double.each, XW)
    multi.obj$int.ma <- exact*multi.obj$int.ma
    multi.obj$int.array <- exact*multi.obj$int.array
    multi.obj$eta <- drop(multi.obj$eta)
    
    int.vec.fix <-  exact * rowSums(multi.obj$help.calc)
    
    Fisher.upp <- matrix(multi.obj$int.array, dimb, dimb)
    
    ## update P
    
    P <- matrix(0,dimb+lin+dimr,dimb+lin+dimr)
    diag(P)[1:nbasis] <- penal * K.base
    
    if(noncat.present)
    {
      alpha.ma.small <- matrix(alpha[(nbasis+1):((m+1)*nbasis)],ncol=m,nrow=nbasis)
      pen.1.vec <- zeta * sqrt(nbasis-1) * apply(alpha.ma.small,2,penal.fct.inv,c.app = c.app, K = K.const) * adaptive.weights[1:m,1]
      pen.2.vec <- (1-zeta) * sqrt(nbasis) * apply(alpha.ma.small,2,penal.fct.inv,c.app = c.app, K = K.smooth) * adaptive.weights[1:m,2]
      for(r in 1:m)
      {
        P[(r*nbasis+1):((r+1)*nbasis),(r*nbasis+1):((r+1)*nbasis)] <-
          xi * ( pen.1.vec[r] * K.const + pen.2.vec[r] * K.smooth) 
      }
    }
    
    if(cat.present)
    {
      levels.ind <- rep(1:number.cat,levels.vec)
      alpha.ma2 <- matrix(alpha[((m+1)*nbasis+1):((m+m2+1)*nbasis)],nbasis,m2)
      eucl.norm.vec <- rep(0,number.cat)
      for(k in 1:number.cat)
        eucl.norm.vec[k] <- sqrt(levels.vec[k]*nbasis)*penal.fct.inv(as.vector(alpha.ma2[,(levels.ind==k)]),K.eucl[[k]],c.app)
      
      pen.1.vec <- zeta * sqrt(nbasis-1) * apply(alpha.ma2,2,penal.fct.inv,c.app = c.app, K = K.const) * adaptive.weights[(m+1):(m+m2),1]
      pen.2.vec <- (1-zeta) * eucl.norm.vec * adaptive.weights[(m+cumsum(levels.vec)),2] 
      
      for(r in 1:m2)
      {
        P[((r+m)*nbasis+1):((r+m+1)*nbasis),((r+m)*nbasis+1):((r+m+1)*nbasis)] <-
          xi * ( pen.1.vec[r] * K.const)
      }
      
      for(r in 1:number.cat)
      {
        P[((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis),((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis)]  <- P[((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis),((cum.levels[r]+m+1)*nbasis+1):((cum.levels[r+1]+m+1)*nbasis)] + xi * ( pen.2.vec[r] * K.eucl[[r]])
      }
      
    }
    
    if(s==1)
    {
      diag(P)[(dimb+lin+1):(dimb+lin+dimr)] <- 1/Q1
    }else{
      Q_inv <- chol2inv(chol(Q1))
      P <- matrix(0,dimb+lin+dimr,dimb+lin+dimr)
      for(jf in 1:n)
        P[(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s),(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s)]<-Q_inv
    }
    
    diag(P)[(nbasis+1):((m+m2+1)*nbasis)] <- diag(P)[(nbasis+1):((m+m2+1)*nbasis)] + ridge.pen * K.ridge
    
    score <- c(delta.Phi - multi.obj$t.eta %*% multi.obj$int.ma,
               t.XW %*%(delta -  (int.vec.fix * multi.obj$eta))) - P %*% Delta[ll,]
    
    Fisher <- matrix(NA,dimb+lin+dimr,dimb+lin+dimr)
    
    Fisher[1:dimb,1:dimb] <- -Fisher.upp
    
    Fisher[(dimb+1):(dimb+lin+dimr),1:dimb] <-  - t.XW %*% (multi.obj$int.ma*multi.obj$eta)
    
    Fisher[1:dimb,(dimb+1):(dimb+lin+dimr)] <- 
      t(Fisher[(dimb+1):(dimb+lin+dimr),1:dimb])
    
    Fisher[(dimb+1):(dimb+lin+dimr),(dimb+1):(dimb+lin+dimr)] <- 
      -t.XW %*% (XW*multi.obj$eta*int.vec.fix)
    
    ###
    InvFisher<-try(chol2inv(chol(P - Fisher)),silent=T)
    if(class(InvFisher)=="try-error")
      InvFisher<-try(solve(P - Fisher),silent=T)
    if(class(InvFisher)=="try-error" || sum(is.na(InvFisher))>0)
      stop("Current Fisher matrix not invertible")  
    
    Delta[ll+1, ] <- Delta[ll, ] + InvFisher %*% score
    
    ### random effects variance
    alpha <- Delta[ll+1,1:dimb]
    if(lin>0)
    {
      fixef <- Delta[ll+1,(dimb+1):(dimb+lin)]
    }else{
      fixef <- NULL
    }
    ranef <- Delta[ll+1,(1+dimb+lin):(dimb+lin+dimr)]
    
    if(s==1)
    {
      Q1 <- sum(diag(InvFisher)[(dimb+1):(dimb+lin+dimr)])+sum(ranef^2)
    }else{ 
      Q1<-InvFisher[(dimb+lin+1):(dimb+lin+s),(dimb+lin+1):(dimb+lin+s)]+ranef[1:s]%*%t(ranef[1:s])
      for (i in 2:n)
        Q1<-Q1+InvFisher[(dimb+lin+(i-1)*s+1):(dimb+lin+i*s),(dimb+lin+(i-1)*s+1):(dimb+lin+i*s)]+ranef[((i-1)*s+1):(i*s)]%*%t(ranef[((i-1)*s+1):(i*s)])
    }  
    Q1<-Q1/n
    
    Q[[ll+1]]<-Q1
    
    ### convergence criteria   
    
    finish<-(sqrt(sum((Delta[ll, ] - Delta[ll+1, ])^2))/sqrt(sum((Delta[ll, ])^2)) < control$conv.eps)
    
    if(finish)
      break
  }
  
  if(s==1)
    Q1<-sqrt(Q1)
  
  
  comp.ma <- matrix(Delta[ll+1,(nbasis+1):((m+m2+1)*nbasis)],nbasis,m+m2)
  which.sel <- (apply(comp.ma,2,eucl.norm))>1
  
  if(xi==0)
  {
    alpha.ma.small <- matrix(alpha[(nbasis+1):((m+1)*nbasis)],ncol=m,nrow=nbasis)
    adaptive.weights <- 1/cbind(apply(alpha.ma.small,2,penal.fct,K=K.const),
                                apply(alpha.ma.small,2,penal.fct,K=K.smooth))
    
    if(cat.present)
    {
      alpha.ma2 <- matrix(alpha[((m+1)*nbasis+1):((m+m2+1)*nbasis)],ncol=m2,nrow=nbasis)
      
      K.eucl <- list()
      for(ee in 1:number.cat)
      {
        K.list <- bdiag(replicate(levels.vec[ee],K.part,simplify = F))
        K.eucl[[ee]] <- as.matrix(crossprod(K.list))
      }
      
      eucl.norm.vec <- rep(0,number.cat)
      for(k in 1:number.cat)
        eucl.norm.vec[k] <-penal.fct(as.vector(alpha.ma2[,(levels.ind==k)]),K.eucl[[k]]) 
      eucl.P <- rep(eucl.norm.vec,levels.vec)
      
      adaptive.weights <- rbind(adaptive.weights,1/cbind(apply(alpha.ma2,2,penal.fct, K = K.const),
                                                         eucl.P))
    }}
  
  time.vary <- as.numeric(alpha[(nbasis+1):((m+m2+1)*nbasis)])
  ## re-transform (centered & standardized) coefficients
  if(control$standardize)
    time.vary <- time.vary/rep(sd.vec,each=nbasis)

  base.final <- as.numeric(center.fct(Delta[ll+1,],mean.vec=mean.vec,sd.vec=sd.vec,nbasis=nbasis,m=m+m2,standardize = control$standardize, center =control$center))
  names(time.vary) <- paste(rep(colnames(U.both)[-1],each=nbasis),rep(1:nbasis,m+m2),sep=".")
  
  colnames(Delta) <-c(paste("baseline",1:nbasis,sep="."),names(time.vary),colnames(X),colnames(W))
  
  coef <- NULL
  if(lin>0)
    coef <- Delta[ll+1,(1+dimb):(dimb+lin)]
  
  ret.obj <- list()
  ret.obj$Delta <- Delta 
  ret.obj$baseline <- base.final
  ret.obj$time.vary <- time.vary
  ret.obj$coefficients <- coef
  ret.obj$ranef <- Delta[ll+1,(1+dimb+lin):(dimb+lin+dimr)]
  ret.obj$nbasis <- nbasis
  ret.obj$spline.degree <- spline.degree
  ret.obj$diff.ord <- diff.ord
  ret.obj$penal <- penal
  ret.obj$Q_long <- Q
  ret.obj$Q <- Q1
  ret.obj$D <- D
  ret.obj$vary.names <- vary.names
  ret.obj$dim.rndeff <- s
  ret.obj$number.effects <- m
  ret.obj$iter <- ll
  ret.obj$time.grid <- time.grid
  ret.obj$exact <- exact
  ret.obj$Phi.big <- Phi.big
  ret.obj$adaptive.weights <- adaptive.weights
  ret.obj$knots <- knots
  ret.obj$rnd <- rnd
  ret.obj$B.unpen.fact <- B.unpen.fact
  ret.obj$Pen <- Pen
  ret.obj$m <- m
  ret.obj$m2 <- m2
  ret.obj$center <- control$center
  ret.obj$standardize <- control$standardize
  ret.obj$mean.vec <- mean.vec
  ret.obj$sd.vec <- sd.vec
  ret.obj$event <- delta
  ret.obj$event.time <- y
  ret.obj$fix <- fix
  ret.obj$vary.coef <- vary.coef
  ret.obj$data <- data
  return(ret.obj)
}

######### Methods

print.pencoxfrail <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nFixed Effects:\n")
  if(!is.null(x$coefficients)){ 
    cat("\nCoefficients:\n") 
    print(x$coefficients) }else{ 
      cat("\nNo time-constant effects included!\n")
    }
  
  cat("\nSmooth (time-varying) Effects:\n")
  if(!is.null(x$time.vary))
  {
    print(x$vary.names)
  }else{
    cat("\nNo time-varying effects included!\n") 
  }  
  
  if(!is.null(x$rnd))
  {  
    cat("\nRandom Effects:\n")
    if(x$dim.rndeff==1) cat("\nStdDev:\n") else cat("\nCov:\n") 
    
    print(x$StdDev)
    
  }else{
    cat("\nNo random effects included!\n")
  }
}

####

summary.pencoxfrail <- function(object, ...)
{
  coef <- object$coefficients
  baseline.coef <- object$baseline
  time.vary <- object$time.vary
  ranef <- object$ranef
  res <- list(call=object$call,
              coefficients=coef,baseline=baseline.coef,
              time.vary=time.vary,ranef=ranef,StdDev=object$StdDev,rnd=object$rnd,
              vary.names=object$vary.names,dim.rndeff=object$dim.rndeff)
  class(res) <- "summary.pencoxfrail"
  res
}

####

print.summary.pencoxfrail <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nFixed Effects:\n")
  if(!is.null(x$coefficients)){ 
    cat("\nCoefficients:\n") 
    print(x$coefficients) }else{ 
      cat("\nNo time-constant effects included!\n")
    }
  
  cat("\nSmooth (time-varying) Effects:\n")
  if(!is.null(x$time.vary))
  {
    print(x$vary.names)
    cat("\nTime-varying effects coefficients:\n")
    print(x$time.vary)
  }else{
    cat("\nNo time-varying effects included!\n") 
  }  
  
  if(!is.null(x$rnd))
  {  
    cat("\nRandom Effects Variance Components:\n")
    if(x$dim.rndeff==1) cat("\nStdDev:\n") else cat("\nCov:\n") 
    print(x$StdDev)
    cat("\nRandom Effects:\n")
    print(x$ranef)
  }else{
    cat("\nNo random effects included!\n")
  }
}

####

plot.pencoxfrail <- function(x,which.comp=NULL,main=NULL,...)
{
  if(!exists("n.grid"))
    n.grid <- 1000
  time.seq <- seq(0,max(x$time.grid),length.out=n.grid)
  
  Design<-bs.design(time.seq, xl = 0, xr = max(time.seq), spline.degree=x$spline.degree, nbasis = x$nbasis)
  Design  <-cbind(Design%*%x$B.unpen.fact,Design%*%x$Pen)
  
  name.vec <- c("baseline",x$vary.names)
  m<-length(name.vec)
  
  smooth.effects <- Design%*%matrix(c(x$baseline,x$time.vary),x$nbasis,m)
  
  ## transform baseline by exp(.)
  smooth.effects[,1] <- exp(smooth.effects[,1])
  
  if(!exists("plot.data"))
    plot.data <- TRUE
  
  if(is.null(main))
    main<-""

  
  if(is.null(which.comp))
    which.comp<-1:m
  
  p<-length(which.comp)
  if(p>9)
    stop("Too many smooth functions! Please specify at maximum nine.")
  
  a<-ceiling(sqrt(p))
  b<-round(sqrt(p))
  if(b==0)
    b<-1

  par(mfrow=c(a,b),mar=c(5, 5.5, 4, 2) + 0.1)
  
  for(i in which.comp)
  {
    plot(time.seq,smooth.effects[,i],type="l", lwd=2, ylab=name.vec[i],
         main=main,xlab="time",...) 
    
     if(plot.data)
        rug(jitter(x$event.time[x$event==1]))
    
  }
}

####

predict.pencoxfrail <- function(object,newdata=NULL,time.grid=NULL,...)
{
  if(is.null(newdata))
    newdata <- object$data 

    if(!is.data.frame(newdata))
    stop("newdata needs to be a data frame!")
  
  data <- newdata  
  
  if(length(object$fix[[2]])==3)
  {
  time.name <- as.character(object$fix[[2]][[2]])
  event.name <- as.character(object$fix[[2]][[3]])
  }else{
  time.name <- as.character(object$fix[[2]][[3]])
  event.name <- as.character(object$fix[[2]][[4]])
  }
  if(is.null(time.grid))
    time.grid <- seq(0,max(data[,time.name]),length.out=1000)
  time.length <- length(time.grid)
  
  data[,time.name][is.na(data[,time.name])] <- max(time.grid)
  data[,event.name][is.na(data[,event.name])] <- 1

  Resp <- model.response(model.frame(object$fix, data))
  
  if(ncol(Resp)==2)
  {
    y0 <- rep(0,nrow(Resp))  
    y <- Resp[,1]
    delta <- Resp[,2]
  }else{
    y0 <- Resp[,1]  
    y <- Resp[,2]
    delta <- Resp[,3]
  }
  N <- length(y)
  
  
  Design<-bs.design(time.grid, xl = 0, xr = max(time.grid), spline.degree=object$spline.degree, nbasis = object$nbasis)
  Design  <-cbind(Design%*%object$B.unpen.fact,Design%*%object$Pen)
  
  smooth.effects <- Design%*%matrix(c(object$baseline,object$time.vary),object$nbasis,object$m+object$m2+1)
  
  vary.names <- attr(terms(object$vary.coef),"term.labels")
  
  ind1 <- sapply(data[,object$vary.names],is.factor)
  ind2 <- sapply(object$vary.names,factor.test)
  ind <- (ind1+ind2)>0
  
  number.cat <- sum(ind)
  number.noncat <- sum(!ind)
  
  cat.present <- number.cat>0
  noncat.present <- number.noncat>0 
  
  if(cat.present)
  {
    vary.coef.cat <- formula(paste("~ 1+",paste(object$vary.names[ind],collapse="+")))
  }else{
    vary.coef.cat <- formula("~ 1")
  }
  
  if(noncat.present)
  {
    vary.coef <- formula(paste("~ 1+",paste(object$vary.names[!ind],collapse="+")))
  }else{
    vary.coef <- formula("~ 1")
  }
  
  if(cat.present)
  {
    levels.vec <- numeric()
    for(i in which(ind==TRUE))
      levels.vec[i] <- nlevels(data[,vary.names[i]])-1
    cum.levels <- c(0,cumsum(levels.vec))
  }
  
  ic.dummy<-attr(terms(object$fix),"intercept")  
  
  if(ic.dummy!=1){
    fix.help <- update(object$fix,~ .+1)
    X <- model.matrix(fix.help, data)[,-1]
  }else{
    X <- model.matrix(object$fix, data)  
  }
  
  ### remove intercept
  X <- X[,-1,drop=F]
  
  ## set up random formula
  rnd.len<-length(object$rnd)
  
  rndformula <- as.character(object$rnd)
  
  trmsrnd <- terms(object$rnd[[1]])
  newrndfrml <- "~ -1"
  newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(object$rnd)[1] else "", sep=" + ")
  
  if(length(attr(trmsrnd, "variables"))>1)
  {
    newrndfrml <- paste(newrndfrml,  
                        paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                          paste(lbl, names(object$rnd)[1], sep=":")
                        }), collapse=" + "), sep="+") 
  }
  
  data[,colnames(data)==(names(object$rnd)[1])] <- as.factor(as.character(data[,colnames(data)==(names(object$rnd)[1])]))
  subject.names<-names(object$rnd)  
  
  if(length(levels(data[,colnames(data)==(names(object$rnd)[1])]))>1)
  {
  W_start <- model.matrix(formula(newrndfrml), data)
  
  rnlabels <- terms(formula(newrndfrml))
  random.labels <- attr(rnlabels,"term.labels")
  id.levels <- levels(data[,colnames(data)==(names(object$rnd)[1])])
  n <- length(id.levels)
  s <- dim(W_start)[2]/n
  
  if(s>1)
  {
    W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
    for (i in 2:n)
      W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
  }else{
    W<-W_start
  }
  
  XW <- cbind(X,W)
  
  ## check if subjects match
  ranef <- rep(0,ncol(W)); names(ranef) <- colnames(W)
  ranef.index <- is.element(names(object$ranef),names(ranef))
  ranef[names(object$ranef[ranef.index])] <- object$ranef[ranef.index]
  ## time-constant part of eta
  eta.vec <- XW %*% c(object$coefficients,ranef)
  }else{
    n <- 1;id.levels <- data[1,colnames(data)==(names(object$rnd)[1])]
    eta.vec <- rep(object$ranef[is.element(names(object$ranef),paste0(names(object$rnd)[1],id.levels))],nrow(data))
    
    if(length(eta.vec)==0)
      eta.vec <- rep(0,nrow(data))
    
    if(!is.null(object$coefficients))
      eta.vec <- eta.vec + X %*% object$coefficients 
  }
  ### time varying coef
  U <- model.matrix(object$vary.coef, data)
  m <- ncol(U)-1
  
  U2 <- model.matrix(vary.coef.cat, data)
  U2 <- U2[,-1]
  m2 <- ncol(U2)
  
  U.both <- cbind(U,U2)
  
  T.mat <- matrix(NA,time.length,N)
  event.ma <- cbind(y0,y)
  T.mat[1:length(time.grid),1:N] <- apply(event.ma,1,function(z,time.grid){(time.grid<z[2] & time.grid>=z[1])},time.grid=time.grid)
  
  
  id.haz.ma <- numeric()
  id.haz.list <- list()
  counter <- 0
  numb.subj.vec <- numeric()
  for(i in 1:n)
  {
    index.help <- (data[,colnames(data)==(names(object$rnd)[1])]==id.levels[i])
    T.mat.help <- T.mat[,index.help,drop=FALSE]
    numb.subj <- sum(T.mat.help[1,])
    U.both.help <- U.both[index.help,,drop=FALSE]
    eta.vec.help <- eta.vec[index.help]
    
    index.vec.help <- which(T.mat.help[1,]==TRUE)
    id.help <- rep(c(1:numb.subj),times=c(diff(index.vec.help),ncol(T.mat.help)-tail(index.vec.help,1)+1))
    
    for(io in 1:numb.subj)
    {  
      id.haz.list[[counter+io]] <- rep(NA,time.length)
      
      U.both.help2 <- U.both.help[(id.help==io),,drop=FALSE]
      eta.vec.help2 <- eta.vec.help[(id.help==io)]
      T.mat.help2 <- T.mat.help[,(id.help==io),drop=FALSE]
      for(j in 1:nrow(U.both.help2))
        id.haz.list[[counter+io]][T.mat.help2[,j]] <- (eta.vec.help2[j]+apply((t(smooth.effects)*U.both.help2[j,]),2,sum))[T.mat.help2[,j]]
  
      id.haz.ma <- cbind(id.haz.ma, exp(id.haz.list[[counter+io]]))
    } 
    counter <- counter + numb.subj
    numb.subj.vec[i] <- numb.subj
  }
  if(all(numb.subj.vec==1))
  {  
  colnames(id.haz.ma) <- paste0(subject.names,id.levels)
  }else{
    label.numb <- sequence(numb.subj.vec)
    colnames(id.haz.ma) <- paste0(rep(paste0(subject.names,id.levels),times=numb.subj.vec),".",label.numb)
  }
  surv.ma <- exp(-(apply(id.haz.ma,2,cumsum)))

 ret.obj <- list()
 ret.obj$time.grid <- time.grid
 ret.obj$haz <- id.haz.ma
 ret.obj$survival <- rbind(1,surv.ma[-nrow(surv.ma),,drop=F])
 return(ret.obj)  
}

