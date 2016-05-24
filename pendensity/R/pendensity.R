pendensity <- function(form,base="bspline",no.base=NULL,max.iter=20,lambda0=500,q=3,sort=TRUE,with.border=NULL,m=q,data=parent.frame(),eps=0.01) {
  #m <- q-1 #order of penalty
  penden.env <- new.env()
  assign("frame",data,penden.env)
  assign("q",q,penden.env)
  assign("m",m,penden.env)
  assign("base",base,penden.env)
  assign("lambda0",lambda0,penden.env)
  assign("form",form,penden.env)
  #accur<- 0.01
  accur <- 0.001
  knots.spline <- NULL
  den.form <-  pendenForm(penden.env)
  stand.num <- NULL
  #eps2 <- 1e-03
  assign("y",y <-den.form$y,penden.env)
  assign("n",n <- length(get("y",penden.env)),penden.env)
  if(sort) {
    y.order <- order(get("y",penden.env))
    assign("y",y<-get("y",penden.env)[y.order],penden.env)
  }
  if(is.null(lambda0)) lambda0 <- 0
  if(!is.null(den.form$no.cov)) {
    assign("N",N <-1,penden.env)
    assign("x.mat",x.mat <- NULL,penden.env)
    assign("x",x <- NULL,penden.env)
    assign("x.null",TRUE,penden.env)
    assign("p",p <- 1,penden.env)
    assign("x.factor",x.factor <- matrix(1,1,1),penden.env)
    assign("len.x.fac",1,penden.env)
    assign("Z",Z <- matrix(1,get("n",penden.env),1),penden.env)
    Z.index <- rep(1,n)
    assign("Z.index",Z.index,penden.env)
  }

  if(!is.null(den.form$parcov$name)) {
    assign("N",N <- 1+length(den.form$parcov$x.mat[1,]),penden.env)
    assign("p",p <- N,penden.env)
    x.1 <- matrix(1,length(den.form$parcov$x.mat[,1]),1)
    assign("x",x <- den.form$parcov$x.mat,penden.env)
    assign("x.null",FALSE,penden.env)
    assign("x.mat",x.mat <- den.form$parcov$x.mat,penden.env)
    assign("Z",Z <- cbind(x.1,den.form$parcov$x.mat),penden.env)
    if(sort) {
      Z <- Z[y.order,]
      assign("Z",Z,penden.env)
    }
    Z.h <- Z
    colnames(Z.h) <- c()
    x.factor <- cbind(matrix(1,length(den.form$parcov$cons[,1]),1),den.form$parcov$cons)
    colnames(x.factor) <- c()
    assign("x.factor",x.factor,penden.env)
    len.x.fac<-length(x.factor[,1])
    assign("len.x.fac",len.x.fac,penden.env)
    Z.index <- c()
    for(i in 1:n) for(j in 1:len.x.fac) if(identical(Z.h[i,],x.factor[j,])) Z.index[i] <- j
    assign("Z.index",Z.index,penden.env)
  }

  if(!is.null(no.base)) {
    if(base=="gaussian") K <- 2*no.base+1
    if(base=="bspline") if(q>=2) K <- 2*no.base+1-q+2
  }
  if(is.null(no.base)) {
    if(base=="gaussian") K <- 41
    if(base=="bspline")  if(q>=2) K <- 41-q+2  else K <- 41
  }
  ymin <- min(get("y",penden.env))
  ymax <- max(get("y",penden.env))

 if(!is.null(with.border)) {
   assign("h",h <- (ymax-ymin)/(K-1-with.border*2),penden.env)
   assign("MeanW",MeanW <- c(seq(ymin-with.border*h,ymin-h,by=h),seq(ymin,ymax,by=h),seq(ymax+h,ymax+with.border*h,by=h)),penden.env)
 }
 else {
   assign("h",h <- (ymax-ymin)/(K-1),penden.env)
   assign("MeanW",MeanW <- seq(ymin,ymax,by=h),penden.env)
 }

  assign("Stand.abw",Stand.abw <- rep((h/2),length(MeanW)),penden.env)
  knots.val<-list(dis="equi",val=MeanW)
  
  nn <- matrix(1:length(y))

  if(base=="gaussian") {
    help.degree <- 0
    base.den <- apply(nn,1,function(i,MeanW,Stand.abw,y) dnorm(y[i],MeanW,Stand.abw),MeanW,Stand.abw,y)
    base.den2 <- NULL
    help.env <- NULL
  }
  
  if(base=="bspline") {
    base.val<- my.bspline(h,q,knots.val,y,K,plot.bsp=FALSE)
    help.env <- base.val$help.env
    assign("base.den",base.den <- base.val$base.den,penden.env)
    base.den2 <- base.val$base.den2
    stand.num <- base.val$stand.num
    assign("knots.val",knots.val <- base.val$knots.val,penden.env)
    assign("K",K <- base.val$K,penden.env)
    help.degree <- base.val$help.degree
  }
  assign("M",M<- floor(K/2) + 1,penden.env)

  liste.opt <- matrix(0,1,5+N*K)
  liste.opt.h <- matrix(0,1,5+N*K)
  if(get("x.null",penden.env)) {
    no.base <- (K-1)/2
    K.help <- seq(-no.base*N,no.base*N,by=1)
    names <- c("lambda0","my.AIC","marg.Likelihood","opt.Likelihood","pen.Likelihood",paste("beta.val.",K.help,sep=""))
    colnames(liste.opt) <- names
  }

  kk <- 1
  listen <- list()
  calc <- TRUE
  L <- L.mat(penden.env)
  D.m(penden.env)
  liste.opt[1,1] <- get("lambda0",penden.env)
  
  while(calc) {
    i <- 1
    liste <- matrix(0,2,2+(N*K))
    n.liste <- matrix(0,1,2+(N*K))
    assign("beta.val",rep(0,K*N),penden.env)
    assign("ck.temp",ck(penden.env,get("beta.val",penden.env)),penden.env)
    assign("f.hat.val",f.hat(penden.env),penden.env)
    f.hat.val <- get("f.hat.val",penden.env)
    liste[i,3:((N*K)+2)] <- get("beta.val",penden.env)
    liste[i,2] <- pen.log.like(penden.env,get("lambda0",penden.env),f.hat.val)      
    i<-2
    val <- new.beta.val(liste[i-1,2],penden.env)
    if(is.na(val$Likelie)) browser()
    if(!is.na(val$Likelie)) {
      liste[i,2] <- val$Likelie
      liste[i,1] <- val$step
      liste[i,3:((N*K)+2)] <- get("beta.val",penden.env)
    }
    if(is.na(val$Likelie) & kk==1) stop ("Conditioning problem. Restart iteration with new start values")
    if(is.na(val$Likelie) & kk!=1) break
    while(liste[i,2]-liste[(i-1),2] > accur) {
      help <- new.beta.val(liste[i,2],penden.env)
      #if(is.na(val$Likelie)) browser()
      if(!is.na(help$Likelie)) {
        i <- i+1
        liste <- rbind(liste,n.liste)
        liste[i,3:((N*K)+2)] <- get("beta.val",penden.env)
        liste[i,2] <- help$Likelie
        liste[i,1] <- help$step
      }
      else break
    }
    help.Derv1 <- Derv1(penden.env)
    help.Derv2 <- Derv2(penden.env,get("lambda0",penden.env))
    assign("Derv1.cal",help.Derv1$Derv1.cal,penden.env)
    assign("f.hat.val",help.Derv1$f.hat.val,penden.env)
    assign("Derv2.pen",help.Derv2$Derv2.pen,penden.env)
    assign("Derv2.cal",help.Derv2$Derv2.cal,penden.env)
    l <- length(liste[,1])
    pen.Likelihood <- liste[l,2]
    opt.Likelihood <-pen.log.like(penden.env,get("lambda0",penden.env))
    marg.likeli <- marg.likelihood(penden.env,pen.Likelihood)
    myAIC.help <- my.AIC(penden.env,lambda0=liste.opt[kk,1],opt.Likelihood)
    my.AIC <- myAIC.help$myAIC
    my.trace <- myAIC.help$mytrace    
    liste.opt[kk,2] <- my.AIC
    liste.opt[kk,6:(N*K+5)] <- get("beta.val",penden.env)
    liste.opt[kk,3] <- marg.likeli
    liste.opt[kk,4] <- opt.Likelihood
    liste.opt[kk,5] <- pen.Likelihood
    if(kk>2) if(abs(liste.opt[kk,3]/liste.opt[kk-1,3]-1)<0.005) {
      calc <- FALSE
      assign("beta.val",liste.opt[kk-1,6:(N*K+5)],penden.env)
      assign("ck.temp",ck(penden.env,get("beta.val",penden.env)),penden.env)
      assign("f.hat.val",f.hat(penden.env),penden.env)
      assign("lambda0",liste.opt[kk-1,1],penden.env)
      #print(get("beta.val",penden.env))
      liste.opt <- liste.opt[-kk,]
      listen <- listen[-kk]
      break
    }
    listen[[kk]] <- liste
    kk <- kk+1
    help.lambda <- new.lambda(penden.env,lambda0=liste.opt[kk-1,1])
    if(abs(help.lambda-liste.opt[kk-1,1])<eps*liste.opt[kk-1,1] | kk>max.iter) calc <- FALSE
    else {
      if(calc) liste.opt <- rbind(liste.opt,liste.opt.h)
      if(calc) assign("lambda0",liste.opt[kk,1] <- help.lambda,penden.env)
    }
  }
 
  lambda0 <- get("lambda0",penden.env)
  opt.Likelihood <-pen.log.like(penden.env,lambda0=0)
  marg.likeli <- marg.likelihood(penden.env,pen.Likelihood)
  beta.val <- matrix(get("beta.val",penden.env),K,N,byrow=TRUE)
  varpar <- variance.par(penden.env)
  biaspar <- bias.par(penden.env)
  
  if(!sort) x.list <- den.form$parcov$x.mat else x.list <- den.form$parcov$x.mat[y.order,]
  
  ck.weights <- ck(penden.env,get("beta.val",penden.env))

  obj <- list(call=form,values=list(y=y,x=x,sort=sort,covariate=list(Z=Z,levels=den.form$parcov$levels,how.levels=den.form$parcov$length.how,how.combi=den.form$parcov$combi.how,x.factor=x.factor)),splines=list(K=K,N=N,MeanW=MeanW,Stand.abw=Stand.abw,h=h,m=m,q=q,stand.num=stand.num,base=base,base.den=base.den,base.den2=base.den2,L=L,Dm=get("Dm",penden.env),help.degree=help.degree,knots.val=knots.val),results=list(ck=ck.weights,beta.val=beta.val,lambda0=lambda0,fitted=get("f.hat.val",penden.env),variance.par=varpar,bias.par=biaspar,AIC=list(my.AIC=my.AIC,my.trace=my.trace),Derv=list(Derv2.pen=get("Derv2.pen",penden.env),Derv2.cal=get("Derv2.cal",penden.env),Derv1.cal=get("Derv1.cal",penden.env),iterations=list(list.opt.results=liste.opt,all.lists=listen)),likelihood=list(pen.Likelihood=pen.Likelihood,opt.Likelihood=opt.Likelihood,marg.Likelihood=marg.likeli)))
  class(obj) <- "pendensity"
  return(obj)
}
