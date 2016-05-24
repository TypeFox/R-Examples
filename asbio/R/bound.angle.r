bound.angle<-function(X,Y,nX,nY){
rX<-nX-X
rY<-nY-Y
theta<-atan(rY/rX)*180/pi
  qrat<-function(rX,rY){
  if(rX==0&rY==0){q=NA} 
  if(rX>0&rY>0){q=1} else
  if(rX<0&rY>0){q=2} else
  if(rX<0&rY<0){q=3} else
  if(rX>0&rY<0){q=4} else
  if(rX==0&rY>0){q=1} else
  if(rX==0&rY<0){q=4} else
  if(rX>0&rY==0){q=3} else
  if(rX<0&rY==0){q=3} 
  q
  }
qm<-matrix(ncol=1,nrow=length(rX))
  for(i in 1: length(rX)){
  qm[i]<-qrat(rX[i],rY[i])
  if(is.na(qm[i])){theta[i]=NA} else
  if(qm[i]==1){theta[i]=theta[i]} else
  if(qm[i]==2){theta[i]=abs(theta[i])+90} else
  if(qm[i]==3&theta[i]!=0){theta[i]=-180+abs(theta[i])} else
  if(qm[i]==3&theta[i]==0){theta[i]=180}
  if(qm[i]==4){theta[i]=theta[i]}
  }
theta<-data.frame(near.angle=theta)
theta
}