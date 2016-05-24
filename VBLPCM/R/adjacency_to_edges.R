Y_to_E<-function(N, NE, directed, Y)
  {
  E<-matrix(NaN,NE,2)
  ans<-.C("Y_to_E", NAOK=TRUE, N=as.integer(N),  
          directed=as.integer(directed), Y=as.numeric(t(Y)), E=as.integer(t(E)),PACKAGE="VBLPCM")
  return(t(matrix(ans$E,2)))
  }
Y_to_nonE<-function(N, NnonE, directed, Y)
  {
  nonE<-matrix(0, NnonE, 2)
  ans<-.C("Y_to_nonE", NAOK=TRUE, N=as.integer(N), directed=as.integer(directed), 
          Y=as.numeric(t(Y)), nonE=as.integer(t(nonE)),PACKAGE="VBLPCM")
  return(t(matrix(ans$nonE,2)))
  }
Y_to_M<-function(N, NM, directed, Y)
  {
  M<-matrix(0, NM, 2)
  ans<-.C("Y_to_M", NAOK=TRUE, N=as.integer(N), directed=as.integer(directed), 
          Y=as.numeric(t(Y)), M=as.integer(t(M)),PACKAGE="VBLPCM")
  return(t(matrix(ans$M,2)))
  }
E_to_Y<-function(N, NE, directed, E)
  {
  Y<-matrix(0,N,N)
  ans<-.C("E_to_Y", NAOK=TRUE, N=as.integer(N), NE=as.integer(NE), 
          directed=as.integer(directed), E=as.integer(t(E)), Y=as.numeric(t(Y)), PACKAGE="VBLPCM")
  return(matrix(ans$Y,N))
  }

hops_to_hopslist<-function(hops,diam,N)
  {
  hopslist<-matrix(0,N,1+diam+N) #nodeid, #hops, nodeids
  for (i in 1:N)
    {
    tmp<-sort(hops[i,],index=1)
    for (h in 0:diam)
      hopslist[i,1+h]<-sum(tmp$x==h)
    hopslist[i,(2+diam):ncol(hopslist)]<-tmp$ix
    }
  return(hopslist)
  }
