clustmember<-function(res,x,mid=T,cl_label="",which=NA,main="Cluster Membership Graph",xlab="Cluster",color=diverge_hcl(101, h = c(0, 130)),...)
{
#require(flexclust)
if(class(res)!="kcca")
{
res<-as.kcca(res,x,...)
}

minx<-min(x)
maxx<-max(x)
mycolor<-color
nx<-dim(res@centers)[1]
ny<-dim(x)[2]
xseq<-seq(0,1,length.out=nx+1)
xticks<-xseq[-length(xseq)]+xseq[2]/2
yseq<-seq(0,1,length.out=ny+1)
yticks<-yseq[-length(yseq)]+yseq[2]/2
midl<-((7*xseq[2])/45)
outl<-(xseq[2]/30)

if(is.na(which)) which<-1:ny

plot.new()
if(!mid)
  {
  for(i in 1:nx)
    {
    for(j in 1:ny)
      {
        identper <- res@cluster==i
        inclustercol<-(round(mean(x[identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        outclustercol<-(round(mean(x[!identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1


        rect(xseq[i]+outl, yseq[j], xticks[i], yseq[j+1], col = mycolor[inclustercol],...)
        rect(xticks[i], yseq[j], xseq[i+1]-outl, yseq[j+1], col = mycolor[outclustercol],...)

      }
    }
  }
else
  {
  for(i in 1:nx)
    {
    for(j in 1:ny)
      {
        identper <- res@cluster==i
        inclustercol<-(round(mean(x[identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        outclustercol<-(round(mean(x[!identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        rect(xseq[i]+outl, yseq[j], xticks[i]-midl, yseq[j+1], col = mycolor[inclustercol],...)
        rect(xticks[i]-midl, yseq[j], xticks[i]+midl, yseq[j+1], col = mycolor[outclustercol],...)
        rect(xticks[i]+midl, yseq[j], xseq[i+1]-outl, yseq[j+1], col = mycolor[inclustercol],...)

      }
    }
  }



#for(i in 1:length(xseq))
#  {
#  lines(c(xseq[i],xseq[i]),c(0,1))
#  }

#for(i in 1:length(yseq))
#  {
#  lines(c(0,1),c(yseq[i],yseq[i]))
#  }

axis(1, at = xticks, labels = paste(cl_label,1:length(xticks)), tick = F,...)

axis(2, at = yticks, labels = colnames(x)[which], tick = F,las=2,...)

axis(4, at = yticks, labels = colnames(x)[which], tick = F,las=2,...)

title(main=main,xlab=xlab,...)
}




######## Biclustermembergraph ####
biclustmember<-function(bicResult,x,mid=T,cl_label="",which=NA,main="BiCluster Membership Graph",xlab="Cluster",color=diverge_hcl(101, h = c(0, 130)),...)
{
minx<-min(x)
maxx<-max(x)
mycolor<-color
nx<-dim(bicResult@NumberxCol)[1]
ny<-dim(bicResult@NumberxCol)[2]
xseq<-seq(0,1,length.out=nx+1)
xticks<-xseq[-length(xseq)]+xseq[2]/2
yseq<-seq(0,1,length.out=ny+1)
yticks<-yseq[-length(yseq)]+yseq[2]/2
midl<-((7*xseq[2])/45)
outl<-(xseq[2]/30)

if(is.na(which)) which<-1:ny

plot.new()
if(!mid)
  {
  for(i in 1:nx)
    {
    for(j in 1:ny)
      {
      if(bicResult@NumberxCol[i,which[j]])
        {
        identper <- bicResult@RowxNumber[,i]
        inclustercol<-(round(mean(x[identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        outclustercol<-(round(mean(x[!identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1


        rect(xseq[i]+outl, yseq[j], xticks[i], yseq[j+1], col = mycolor[inclustercol],...)
        rect(xticks[i], yseq[j], xseq[i+1]-outl, yseq[j+1], col = mycolor[outclustercol],...)
        }
      }
    }
  }
else
  {
  for(i in 1:nx)
    {
    for(j in 1:ny)
      {
      if(bicResult@NumberxCol[i,which[j]])
        {
        identper <- bicResult@RowxNumber[,i]
        inclustercol<-(round(mean(x[identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        outclustercol<-(round(mean(x[!identper,which[j]]),2)-minx)*(100/(maxx-minx)) +1
        rect(xseq[i]+outl, yseq[j], xticks[i]-midl, yseq[j+1], col = mycolor[inclustercol],...)
        rect(xticks[i]-midl, yseq[j], xticks[i]+midl, yseq[j+1], col = mycolor[outclustercol],...)
        rect(xticks[i]+midl, yseq[j], xseq[i+1]-outl, yseq[j+1], col = mycolor[inclustercol],...)
        }
      else
        {
        rect(xseq[i]+outl, yseq[j], xseq[i+1]-outl, yseq[j+1])
        }
      }
    }
  }



#for(i in 1:length(xseq))
#  {
#  lines(c(xseq[i],xseq[i]),c(0,1))
#  }

#for(i in 1:length(yseq))
#  {
#  lines(c(0,1),c(yseq[i],yseq[i]))
#  }

axis(1, at = xticks, labels = paste(cl_label,1:length(xticks)), tick = F,...)

axis(2, at = yticks, labels = colnames(x)[which], tick = F,las=2,...)

axis(4, at = yticks, labels = colnames(x)[which], tick = F,las=2,...)

title(main=main,xlab=xlab,...)
}


