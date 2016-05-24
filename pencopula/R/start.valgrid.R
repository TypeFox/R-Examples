start.valgrid <- function(penden.env) {
  
  p <- get("p",penden.env)
  X.knots <- matrix(NA,get("DD",penden.env),p)
  tilde.Psi.d.knots.start.r <-  array(NA, dim=c(dim(X.knots)[1],get("ddb",penden.env),p))
    
  for(j in 1:p)  X.knots[,j]  <- get("knots.t",penden.env)[get("Index.basis.D",penden.env)[,j]]
  env.extend <- list()
  length.cond <- get("ddb",penden.env)
  
  for(j in 1:p) {
    name <- c(paste("y",j,sep=""))
    env.extend[[noquote(name)]] <- seq(0,1,length=length.cond)
  }
  assign("X.knots.g.all",expand.grid(env.extend),penden.env)
  
  if(get("adapt.grid",penden.env)) grid.points(penden.env)
  
  if(get("adapt.grid",penden.env)) tilde.Psi.d.knots.start.g <-  array(NA, dim=c(dim(get("X.knots.g",penden.env))[1],get("ddb",penden.env),p))
  tilde.Psi.d.knots.start.g.all <- array(NA, dim=c(dim(get("X.knots.g.all",penden.env))[1],get("ddb",penden.env),p))
  
  assign("X.knots",X.knots,penden.env)

  for (j in 1:p)
    {
      tilde.Psi.d.knots.start.r[,,j] <-  hierarch.bs(X.knots[,j], d = get("d",penden.env), plot.bsp = FALSE,typ=3,penden.env,int=FALSE)$B.tilde
      if(get("adapt.grid",penden.env)) tilde.Psi.d.knots.start.g[,,j] <-  hierarch.bs(get("X.knots.g",penden.env)[,j], d = get("d",penden.env), plot.bsp =FALSE,typ=3,penden.env,int=FALSE)$B.tilde
      tilde.Psi.d.knots.start.g.all[,,j] <-  hierarch.bs(get("X.knots.g.all",penden.env)[,j], d = get("d",penden.env), plot.bsp =FALSE,typ=3,penden.env,int=FALSE)$B.tilde
    }

  assign("tilde.PSI.d.D.knots.start.r",tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,1],1],penden.env)
  if(get("adapt.grid",penden.env)) assign("tilde.PSI.d.D.knots.start.g",tilde.Psi.d.knots.start.g[,get("Index.basis.D",penden.env)[,1],1],penden.env)
  assign("tilde.PSI.d.D.knots.start.g.all",tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,1],1],penden.env)

  
  for (j in 2:p)
    {
      assign("tilde.PSI.d.D.knots.start.r",get("tilde.PSI.d.D.knots.start.r",penden.env) * tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,j],j],penden.env)
      if(get("adapt.grid",penden.env)) assign("tilde.PSI.d.D.knots.start.g",get("tilde.PSI.d.D.knots.start.g",penden.env) * tilde.Psi.d.knots.start.g[,get("Index.basis.D",penden.env)[,j],j],penden.env)
      assign("tilde.PSI.d.D.knots.start.g.all",get("tilde.PSI.d.D.knots.start.g.all",penden.env) * tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,j],j],penden.env)
    }
  if(get("base",penden.env)=="B-spline") assign("ck.val",matrix(solve(get("tilde.PSI.d.D.knots.start.r",penden.env),rep(1,get("DD",penden.env)))),penden.env)
  if(get("base",penden.env)=="Bernstein") assign("ck.val",matrix(rep(1/get("DD",penden.env)),get("DD",penden.env)),penden.env)
  #if((get("base",penden.env)=="Bernstein") & (p>2)) assign("ck.val",matrix(solve(get("tilde.PSI.d.D.knots.start.r",penden.env),rep(1,get("DD",penden.env)))),penden.env)
}
