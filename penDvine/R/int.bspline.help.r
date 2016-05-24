int.bspline.help <- function(penden.env) {
  len.k <- get("K",penden.env)
  x.help <- seq(0,1,length=501)
  int.base <- distr.func.help(base=get("base",penden.env),knots=seq(0,1,length=(get("K",penden.env)+get("q",penden.env)-1)),penden.env,q=get("q",penden.env),y=seq(0,1,length=501))
  assign("int.base",int.base,penden.env)
  assign("x.help",x.help,penden.env)
}
