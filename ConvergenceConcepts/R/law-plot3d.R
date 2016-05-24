law.plot3d <- function(data,probfunc,tinf=-5,tsup=5) {

  nmax <- ncol(data)

  Fhat<-function(x) {
    n<-x[1]
    t<-x[2]
    Fn<-ecdf(data[,n])
    Fn(t)
  }

  g <- expand.grid(t = seq(from=tinf,to=tsup,length=10), n = 1:nmax)
  
  Mymat <- g[,2:1]

  g$z <- abs(apply(Mymat,FUN=Fhat,MARGIN=1)-probfunc(g$t))
  
  print(wireframe(z ~ t + n, data = g,scales = list(arrows = FALSE),drape = TRUE, colorkey = TRUE,screen = list(z = -50, x = -70),zlab=list(expression(hat(l)[n]~'(t)=|'~hat(F)[n]~'(t)-'~F~'(t)|'),rot=90) ,main="Convergence in law?"))
  

}
