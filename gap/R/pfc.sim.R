pfc.sim <- function(famdata,n.sim=1000000,n.loop=1)
{
  famsize<-dim(famdata)[1]
  obsp<-0
  tailpl<-tailpu<-rep(0,n.loop)
  z<-.Fortran("runifamily",famdata=as.integer(matrix(famdata,ncol=3)),famsize=as.integer(famsize),
               nsim=as.integer(n.sim),ncycle=as.integer(n.loop),
               obsp=as.double(obsp),tailpl=as.double(tailpl),tailpu=as.integer(tailpu),PACKAGE="gap")

  list(n.sim=n.sim,n.loop=n.loop,p=z$obsp,tailpl=z$tailpl,tailpu=z$tailpu/n.sim)

}
