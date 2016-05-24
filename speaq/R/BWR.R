BWR <-function(X, groupLabel)
{
  groupNum=length(levels(groupLabel));
  m=colMeans(X);
  groupMean=list();
  groupW=list();
  B=X;
  B[,]=0;
  for (i in 1:groupNum){
    groupLabeli=which(groupLabel==levels(groupLabel)[i]);
    Xi=X[groupLabeli,]
    mi=colMeans(Xi);
    groupMean[[i]]=mi;
    tempi=matrix(data=rep(mi,nrow(Xi)),nrow=nrow(Xi),
        ncol=length(mi),byrow=TRUE);

    Wi=(Xi-tempi)^2; 
    groupW[[i]]=Wi;
    temp=matrix(data=rep(m,nrow(Xi)),nrow=nrow(Xi),
        ncol=length(m),byrow=TRUE);
    B[groupLabeli,]=(tempi-temp)^2;    
  }  
  
  BW=double(ncol(B));
  for (i in 1:length(BW)){
    BW_denominator=0;
    for (j in 1:groupNum){
      BW_denominator=BW_denominator+sum(groupW[[j]][,i]);
    }
    Bw_numerator=sum(B[,i]);
    BW[i]=Bw_numerator/BW_denominator;    
  }
  return(BW); 
}