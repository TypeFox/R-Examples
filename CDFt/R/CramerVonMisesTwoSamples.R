CramerVonMisesTwoSamples <- function(S1, S2){

  xS1=sort(S1)
  M=length(xS1)
  xS2=sort(S2)
  N=length(xS2)

  a=data.frame(val=xS1,rang=seq(M),ens=rep(1,M))
  b=data.frame(val=xS2,rang=seq(N),ens=rep(2,N))
  c=rbind(a,b)


  c=c[order(c$val),]
  c=data.frame(c,rangTot=seq(M+N))

  dtfM=c[which(c$ens==1),]
  dtfN=c[which(c$ens==2),]

  somN = sum( (dtfN$rang - dtfN$rangTot)**2 )
  somM = sum( (dtfM$rang - dtfM$rangTot)**2 )

  U = N*somN + M*somM

  #CvM = (U / (N*M*(N+M))) - ((4*M*N - 1)/(6*(M+N)))
  CvM = (  (U / (N*M)) / (N+M) ) - ((4*M*N - 1)/(6*(M+N)))

  return(CvM)

}
