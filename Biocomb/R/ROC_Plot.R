CalculateHUM_Plot<-function(sel,Sn,Sp,optSn,optSp,HUM,print.optim=TRUE)
{
  x_coor=1-Sn[[sel]]
plot(x_coor, Sp[[sel]], xlab="1 - Specificity", ylab="Sensitivity", type="n", xlim=c(0,1), ylim=c(0,1), lwd=2, asp=1)
polygon(c(1,x_coor, 0,1), c(0,Sp[[sel]], 0,0), col="gainsboro",angle=45,border=NULL)
text(0.5,0.5, sprintf("AUC= %.3f", HUM[sel]), adj=c(0,1), cex=1, col="red")
if(print.optim)
{
points(1-optSn[sel],optSp[sel],cex=1, col="blue")
text(1-optSn[sel],optSp[sel], sprintf("(%.3f,%.3f)",optSp[sel],optSn[sel]), adj=c(0,1), cex=1, col="blue")
}
lines(x_coor, Sp[[sel]], type="l", lwd=2, col="red")
}

Calculate3D<-function(sel,Sn,Sp,S3,optSn,optSp,optS3,thresholds,HUM,name,print.optim=TRUE)
{
  len=length(thresholds)
  z=matrix(S3[[sel]],len,len,byrow=TRUE)
  ym=matrix(Sp[[sel]],len,len,byrow=TRUE)

  vrem=seq(1,length(Sn[[sel]]),by=len)
  x=Sn[[sel]][vrem]
  y=Sp[[sel]][1:len]
  xrow=unique(x)
  yrow=unique(y)

  zz=matrix(0,length(xrow),length(yrow))

  indexrow=NULL
  indexcol=NULL
  for(i in 1:length(xrow))
  {
    rr=which(x==xrow[i])
    zr=z[rr,]
    for(j in 1:length(yrow))
    {
      index=which(ym[rr,]==yrow[j])
      if(length(index)>0)
      {
      zz[i,j]=max(zr[index])
      }
    }
  }
  for(j in 1:length(xrow))
  {
  # insert to exclude zeros in matrix zz
  index=which(zz[j,]==0)
  k=1
  if(length(index)!=0)
  {
  if(length(index)!=1)
  {
    for(i in 2:length(index))
    {
      if(index[i]!=(index[i-1]+1))
      {
        for(ll in k:(i-1))
        {
          zz[j,index[ll]]=zz[j,index[i-1]+1]
        }
        k=i
      }
    }
    if(index[length(index)]<length(yrow))
    {
    for(ll in k:length(index))
    {
      zz[j,index[ll]]=zz[j,index[length(index)]+1]
    }
    }
  }
  else
  {
    if(index[1]<length(yrow))
    {
      zz[j,index[1]]=zz[j,index[1]+1]
    }
  }
  }
  }
  #--------
  if(length(which(loadedNamespaces()=="rgl"))!=0)
  {
    clear3d()
  out=persp3d(xrow,yrow,zz,theta = 120, phi = 10, expand = 0.5,ticktype = "detailed",col="#CC00FFFF",
          ltheta = -120, shade = 0.75, border = NA, main=sel,xlab = name[1], ylab = name[2], zlab = name[3],
          xlim=c(0,1), ylim=c(0,1), zlim=c(0,1))
  if(print.optim)
  {
  points3d(optSn[sel],optSp[sel],optS3[sel],col="red",size=8)
  text3d(optSn[sel],optSp[sel],optS3[sel]+0.2,sprintf("(%.3f,%.3f,%.3f)",optSn[sel],optSp[sel],optS3[sel]),col="blue")
  }
  rgl.bringtotop()
  }
  persp(xrow,yrow,zz,theta = 50, phi = 15, expand = 0.5,ticktype = "detailed",col="#CC00FFFF",
        ltheta = -120, shade = 0.75, border = NA, main=sel,xlab = name[1], ylab =name[2], zlab = name[3],xlim=c(0,1), ylim=c(0,1), zlim=c(0,1),
        cex.lab=2, cex.main=2, cex.axis=1.5)

}

