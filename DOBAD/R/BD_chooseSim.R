##Chooses between different methods of simulation.
## use the ap
chooseSim.condBD.1 <- function(T = 1, a = 5, b = 7, L = 0.5, m = 0.7, nu = 0.4, n.fft = 1024,
                               maxARsims=100) {
  eps=1/maxARsims; ##
  p <- process.prob.one(T,L,m,nu,X0=a,Xt=b)
  if (p < eps){
    ##print("A")
    res <- sim.condBD.1(T,a,b,L,m,nu,n.fft)
  }
  else{
    ##print("B")
    res <- bdARsimCondEnd.1(Naccepted=1, Nmax=maxARsims, T=T,a=a,b=b,L=L,m=m,nu=nu)
    if(length(res)==0){ ## didn't get a sim out of it.
      ##print("C")
      res <- sim.condBD.1(T,a,b,L,m,nu,n.fft)
    }
    else
      res <- res[[1]]
  }
  return(res)
}
