fruchterman_reingold<-function(net, D=2, steps=1000, repulserad=N^D, m=N*(D-1), volume=N*(D-1))
  {
  directed<-is.directed(net)
  N<-network.size(net)
  Y<-as.sociomatrix(net)
  X<-matrix(runif(N*D,-D^2,D^2),ncol=D)
  C_fruchterman_reingold<- function(directed, N, D, steps, Y, X, repulserad, m, volume)
     {
     ans<-.C("fruchterman_reingold", NAOK=TRUE, directed=as.integer(directed), N=as.integer(N), 
             D=as.integer(D), steps=as.integer(steps), Y=as.double(t(Y)), X=as.numeric(t(X)), 
	     repulserad=as.double(repulserad), m=as.numeric(m), volume=as.numeric(volume),
	     PACKAGE="VBLPCM")
             return(ans)
     }
  
  # update X
  out<-C_fruchterman_reingold(directed, N, D, steps, Y, X, repulserad, m, volume)
  X<-t(matrix(out$X,ncol=N))
  # centre
  X <- X - t(matrix(rep(apply(X,2,mean),N),nrow=D))
  # scale to lie in -5 to 5
  #X<-t(t(X)/apply(X,2,sd)*5)
  return (X)
  }
