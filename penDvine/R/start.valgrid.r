start.valgrid <- function(penden.env) {
  p <- get("p",penden.env)
  X.knots <- matrix(NA,get("DD",penden.env),p)
  tilde.Psi.d.knots.start.r <-  array(NA, dim=c(dim(X.knots)[1],get("ddb",penden.env),p))
  for(j in 1:p)  X.knots[,j]  <- get("knots",penden.env)[get("Index.basis.D",penden.env)[,j]]
  env.extend <- list()
  length.cond <- get("ddb",penden.env)
  for(j in 1:p) {
    name <- c(paste("y",j,sep=""))
    env.extend[[noquote(name)]] <- seq(0,1,length=length.cond)
  }
  assign("X.knots.g.all",expand.grid(env.extend),penden.env)
  tilde.Psi.d.knots.start.g.all <- array(NA, dim=c(dim(get("X.knots.g.all",penden.env))[1],get("ddb",penden.env),p))
  assign("X.knots",X.knots,penden.env)
  index.b <- matrix(0:(get("ddb",penden.env)-1))
  for (j in 1:p)
    {
      if(get("base",penden.env)=="Bernstein") {
        tilde.Psi.d.knots.start.r[,,j] <-  apply(index.b,1,bernstein,x=X.knots[,j],n=get("dd",penden.env))
        tilde.Psi.d.knots.start.g.all[,,j] <-  apply(index.b,1,bernstein,x=get("X.knots.g.all",penden.env)[,j],n=get("dd",penden.env))
      }
      if(get("base",penden.env)=="B-spline") {
        tilde.Psi.d.knots.start.r[,,j] <-  my.bspline(y=X.knots[,j],K=get("K",penden.env)+get("q",penden.env)-1,q=get("q",penden.env))$base.den
        tilde.Psi.d.knots.start.g.all[,,j] <-  my.bspline(y=get("X.knots.g.all",penden.env)[,j],K=get("K",penden.env)+get("q",penden.env)-1,q=get("q",penden.env))$base.den
      }
    }
  assign("tilde.PSI.d.D.knots.start.r",tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,1],1],penden.env)
  assign("tilde.PSI.d.D.knots.start.g.all",tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,1],1],penden.env)
  for (j in 2:p) {
    assign("tilde.PSI.d.D.knots.start.r",get("tilde.PSI.d.D.knots.start.r",penden.env) * tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,j],j],penden.env)
    assign("tilde.PSI.d.D.knots.start.g.all",get("tilde.PSI.d.D.knots.start.g.all",penden.env) * tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,j],j],penden.env)
   }
  
  assign("ck.val",matrix(solve(get("tilde.PSI.d.D.knots.start.r",penden.env),rep(1,get("DD",penden.env)),tol=1e-50)),penden.env)
  assign(x="ck.val.temp",value=get("ck.val",penden.env),envir=penden.env)
}
