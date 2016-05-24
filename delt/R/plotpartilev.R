plotpartilev<-function(pa,dendat=NULL,restri=NULL,pch=21,col="blue")
{
recs<-pa$recs

if (!is.null(dendat)) plot(dendat,xlab="",ylab="",pch=pch)
else{
  xmin<-min(recs[,1])
  xmax<-max(recs[,2])
  ymin<-min(recs[,3])
  ymax<-max(recs[,4])
  plot(0,0,type="n",ylim=c(ymin,ymax),xlab="",ylab="",xlim=c(xmin,xmax))
}

if (is.null(dim(recs))) len<-1 else len<-dim(recs)[1]

if (is.null(restri)) restric<-matrix(1,len,1)
else{
  restric<-matrix(0,len,1)
  restrilen<-length(restri)
  for (i in 1:restrilen){
     restric[restri[i]]<-1
  }
}

i<-1
while (i<=len){

 if (restric[i]==1){

    if (len==1){
      x<-c(recs[1],recs[1],recs[2],recs[2])
      y<-c(recs[3],recs[4],recs[4],recs[3])
    }
    else{
      x<-c(recs[i,1],recs[i,1],recs[i,2],recs[i,2])
      y<-c(recs[i,3],recs[i,4],recs[i,4],recs[i,3])
    }

    if (pa$values[i]==0) colo<-NA else colo=col  
    polygon(x,y,col=colo)

    #lines(c(recs[i,1],recs[i,1]),c(recs[i,3],recs[i,4]))
    #lines(c(recs[i,1],recs[i,2]),c(recs[i,4],recs[i,4]))
    #lines(c(recs[i,2],recs[i,2]),c(recs[i,4],recs[i,3]))
    #lines(c(recs[i,2],recs[i,1]),c(recs[i,3],recs[i,3]))
 }

 i<-i+1
}

if (!is.null(dendat)){

  points(dendat,pch=pch)

  suppor<-supp(dendat)

  x<-c(suppor[1],suppor[1],suppor[2],suppor[2])
  y<-c(suppor[3],suppor[4],suppor[4],suppor[3])
  polygon(x,y)
  #lines(c(suppor[1],suppor[1]),c(suppor[3],suppor[4]))
  #lines(c(suppor[1],suppor[2]),c(suppor[4],suppor[4]))
  #lines(c(suppor[2],suppor[2]),c(suppor[4],suppor[3]))
  #lines(c(suppor[2],suppor[1]),c(suppor[3],suppor[3]))
}

}






