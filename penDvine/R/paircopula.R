paircopula <- function(data,K=8,base="Bernstein",max.iter=30,lambda=100,data.frame=parent.frame(),m=2,
                       fix.lambda=FALSE,pen=1,q=2) {
  penden.env <- new.env()
  assign("m",m,penden.env)
  assign("d",K,penden.env)
  assign("n",dim(data)[1],penden.env)
  assign("frame",data.frame,penden.env)
  assign("Y",data,penden.env)
  assign("D",D<-K,penden.env)
  assign("no",FALSE,penden.env)
  assign("max.iter",max.iter,penden.env)
  assign("base",base,penden.env)
  dd <- K # Grad der Bernsteinbasis
  assign("dd",dd,penden.env)
  if(base=="Bernstein") ddb <- dd+1 #Anzahl Bernsteinpolynome
  if(base=="B-spline" & q==2) ddb <- K+q-1
  if(base=="B-spline" & q==1) ddb <- K+q-2
  assign("p",p <- 2,penden.env)
  assign("fix.lambda",fix.lambda,penden.env)
  assign("pen",pen,penden.env)
  assign("q",q,penden.env)
  assign("base",base,penden.env)
  assign("lambda",lambda,penden.env)
  DD <- ddb^2
  assign("knots",seq(0,1,length=ddb),penden.env)
  
  Index.basis.D <- matrix(NA,DD,2)
  Index.basis.D[,1] <- rep(seq(1,ddb),ddb)
  Index.basis.D[,2] <- sort(Index.basis.D[,1])
  assign("Index.basis.D",Index.basis.D,penden.env)
  
  T.marg <- array(NA, dim=c(ddb,DD,p))

    for ( j in 1:p)
    {
      for ( l in 1:ddb)
        {
          T.marg[l,,j] <- (Index.basis.D[,j] == l)+0
        }
    }
  
  assign("T.marg",T.marg,penden.env) 
  tilde.Psi.d <- array(NA,dim=c(get("n",penden.env),ddb,p))
  index.b <- matrix(0:dd)

  if(base=="Bernstein") for(j in 1:p) tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=get("Y",penden.env)[,j],n=dd)
  if(base=="B-spline") for(j in 1:p) tilde.Psi.d[,,j] <- my.bspline(y=get("Y",penden.env)[,j],K=K+q-1,q=q)$base.den
  if(base=="B-spline") {
     help <- my.bspline(y=seq(0,1,length=ddb),K=K+q-1,q=q)
     assign("C",help$base.den,penden.env)
     assign("stand.num",help$stand.num,penden.env)
  }
  
  assign("tilde.Psi.d",tilde.Psi.d,penden.env)

  assign("tilde.PSI.d.D",tilde.Psi.d[,Index.basis.D[,1],1],penden.env)

  for (j in 2:p)
    {
      assign("tilde.PSI.d.D",get("tilde.PSI.d.D",penden.env) * get("tilde.Psi.d",penden.env)[,Index.basis.D[,j],j],penden.env)
    }

  start.valgrid(penden.env)
  
  if(base=="Bernstein") int.bernstein.help(penden.env)
  if(base=="Bernstein") assign("tilde.Psi.knots.d",apply(index.b,1,bernstein,x=get("knots",penden.env),n=dd),penden.env)
  if(base=="B-spline") assign("tilde.Psi.knots.d", my.bspline(y=get("knots",penden.env),K=K+q-1,q=q)$base.den,penden.env)
  if(base=="B-spline") int.bspline.help(penden.env)
    
  A <- array(NA, dim=c(get("ddb",penden.env),DD,p))
  
  for ( j in 1:p)
    {
      A[,,j] <- get("tilde.Psi.knots.d",penden.env) %*% T.marg[,,j]
    }

  assign("A.Restrict",A,penden.env)
  pen.matrix(penden.env)

  liste <- matrix(0,1,3+DD+1+1)
  n.liste <- matrix(0,1,3+DD+1+1)
  lam <- coef <- c()
  lam <- "lambda1"
  for(j in 1:DD) coef[j] <- paste("b.",j,sep="")
 
  colnames(liste) <- c("pen.log.like","log.like","marg.log.like",lam,"cAIC",coef)
  help.str <- paste("d=",get("d",penden.env),"D=",get("D",penden.env),"lambda=",get("lambda",penden.env)[1],sep="")
  assign("help.str",help.str,penden.env)

  f.hat.val(penden.env,cal=TRUE)
  if(get("no",penden.env)) {
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("cAIC",0,penden.env)
    assign("BIC",0,penden.env)
    class(penden.env) <- "pencopula"
    return(penden.env)
  }
  pen.log.like(penden.env,cal=TRUE)
  Derv1(penden.env)
  Derv2(penden.env)
  marg.likelihood(penden.env,pen.likelihood=get("pen.log.like",penden.env))
  my.IC(penden.env)

  assign("i",i <- 1,penden.env)
  liste[i,1] <- get("pen.log.like",penden.env)
  liste[i,2] <- get("log.like",penden.env)
  liste[i,3] <- get("marg.log.like",penden.env)
  liste[i,4] <- get("lambda",penden.env)
  liste[i,5] <- get("cAIC",penden.env)
  liste[i,(6:(6+DD-1))] <- get("ck.val",penden.env)

  assign("liste",liste,penden.env)
  f.hat.val(penden.env,temp=TRUE)
  if(new.weights(penden.env,lambda.temp=lambda)=="fehler"){
    assign("pen.log.like",0,penden.env)
    assign("log.like",0,penden.env)
    assign("cAIC",0,penden.env)
    assign("BIC",0,penden.env)
    class(penden.env) <- "paircopula"
    return(penden.env)
  }
  if(get("lambda.change",penden.env)) {
    pen.log.like(penden.env,cal=TRUE)
    Derv1(penden.env)
    Derv2(penden.env)
    marg.likelihood(penden.env,pen.likelihood=get("pen.log.like",penden.env))
    my.IC(penden.env)
    liste[i,1] <- get("pen.log.like",penden.env)
    liste[i,2] <- get("log.like",penden.env)
    liste[i,3] <- get("marg.log.like",penden.env)
    liste[i,4] <- get("lambda",penden.env)
    liste[i,5] <- get("cAIC",penden.env)
    assign("liste",liste,penden.env)
  }

  my.loop(penden.env)

  my.IC(penden.env)
  class(penden.env) <- "paircopula"
  return(penden.env)
}
 

