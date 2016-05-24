# 'modify Z in sparse format for position (i,j)
modsparse_cpp<-function(Zsparse=Zsparse,i=i,j=j,relax=F,rmax=NULL){
  if(is.null(rmax)){
    rmax=Zsparse$p
  }
  res=.Call( "modsparse_cpp",Zsparse$Zi,Zsparse$Zj,Zsparse$Si,Zsparse$Sj,Zsparse$compl,Zsparse$p,i,j,relax,rmax, PACKAGE = "CorReg")
  
  if(res$mod){
    Zsparse$Zi=res$Zi
    Zsparse$Zj=res$Zj
    Zsparse$Si=res$Si
    Zsparse$Sj=res$Sj
    Zsparse$compl=res$complexite
  }
  
  return(list(Zsparse=Zsparse,mod=res$mod))
  
}