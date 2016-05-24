cond.cop <- function(data,coef,K,diff="u2",Index.basis.D,base,q=2) {
  p <- 2
  penden.env <- new.env()
  assign("K",K,penden.env)
  assign("base",base,penden.env)
  assign("q",q,penden.env)
  if(base=="Bernstein") int.bernstein.help(penden.env)
  if(base=="B-spline") int.bspline.help(penden.env)
  ddb <- K+1
  if(base=="B-spline"&q==2) ddb <- K-1+q
  if(base=="B-spline"&q==1) ddb <- K-2+q
  assign("u1",data[,1],penden.env)
  assign("u2",data[,2],penden.env)
  tilde.Psi.d <-  array(NA, dim=c(dim(data)[1],ddb,p))
  index.b <- matrix(0:K)
  for (j in 1:p)
    {
      obj <- paste("u",j,sep="")
      if(base=="Bernstein") {
        if(obj==diff) tilde.Psi.d[,,j] <-  apply(index.b,1,bernstein,get(paste("u",j,sep=""),penden.env), n=K)
        if(obj!=diff) tilde.Psi.d[,,j] <-  int.bernstein(penden.env,Y=get(paste("u",j,sep=""),penden.env))
      }
      if(base=="B-spline") {
        if(obj==diff) tilde.Psi.d[,,j] <-  my.bspline(y=get(paste("u",j,sep=""),penden.env),K=K+q-1,q=q)$base.den
        if(obj!=diff) tilde.Psi.d[,,j] <-  int.bspline2(penden.env,Y=get(paste("u",j,sep=""),penden.env))
      }
    }
  assign("tilde.Psi.d",tilde.Psi.d,penden.env)
  assign("tilde.PSI.d.D",tilde.Psi.d[,Index.basis.D[,1],1],penden.env)
  
  for (j in 2:p)
    {
      assign("tilde.PSI.d.D",get("tilde.PSI.d.D",penden.env) * get("tilde.Psi.d",penden.env)[,Index.basis.D[,j],j],penden.env)
    }

  return(get("tilde.PSI.d.D",penden.env)%*%coef)
}
