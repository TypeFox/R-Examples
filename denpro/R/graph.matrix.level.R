graph.matrix.level<-function(dendat,tt=NULL,permu=seq(1:dim(dendat)[1]),
col=seq(1:2000),config="new",shift=0.1,segme=TRUE,poin=FALSE,epsi=0,
ystart=0.5,pch=21,cex=1, yaxt="s",cex.axis=1,texto=TRUE)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

origins<-matrix(0,d,1)
starts<-matrix(0,d,1)
range<-matrix(0,d,1)
minis<-matrix(0,d,1)
means<-matrix(0,d,1)
for (i in 1:d){
    minis[i]<-min(dendat[,i]) 
    means[i]<-mean(dendat[,i])
}

for (i in 1:d) range[i]<-max(dendat[,i])-min(dendat[,i]) 
starts[1]<-0 #min(dendat[,1])
i<-2
while (i<=d){ 
      starts[i]<-starts[i-1]+range[i-1]+epsi
      i<-i+1
}
if (config=="new")
  for (i in 1:d) origins[i]<-starts[i]+mean(dendat[,i])-min(dendat[,i])
else
  for (i in 1:d) origins[i]<-starts[i]+range[i]/2-min(dendat[,i])

  #starts[i]+range[i]/2

plot(x="",y="",xlim=c(starts[1],starts[d]+range[d]),ylim=c(ystart,n+0.5),
xlab="",ylab="",xaxt="n",yaxt=yaxt,cex.axis=cex.axis)

if (is.null(tt)){

   for (j in 1:d){
      for (i in 1:n){
           indo<-permu[i]
           x0<-dendat[indo,j]+starts[j]-min(dendat[,j]) 
                    #dendat[indo,j]+origins[j]
           y0<-i
           x1<-origins[j]  
           y1<-i
           if (segme)
           segments(x0, y0, x1, y1, col=col[indo])
           if (poin)          
           points(x0, y0, col=col[indo], pch=pch, cex=cex)


       }
       if (j>1){
            beg<-starts[j]-epsi/2 
            segments(beg,0.5,beg,n+0.5,lty=1,lwd=2)
       }
       segments(origins[j],1,origins[j],n,lty=2)
       if (texto){
       text(origins[j],ystart,as.character(round(means[j],digits=2)),cex=cex)
       text(starts[j]+shift,ystart,as.character(round(minis[j],digits=2)),
       cex=cex)
       }
   }

}

else{
   paletti<-seq(1:2000)
   coli<-colobary(tt$parent,paletti)  #,segtype="num")

   cente<-c(mean(dendat[,1]),mean(dendat[,2]))
   dendat3<-matrix(0,n,d)
   newcolo<-matrix(0,n,1)
   maxseg<-max(coli)
   curbeg<-1
   i<-maxseg
   while (i >= 1){
       inditree<-which(coli==i)
       indidendat<-tt$infopointer[inditree]
       curend<-curbeg+length(indidendat)-1
       curseg<-dendat[indidendat,]
       leveli<-sqrt(sum(curseg-cente)^2)    #pituus(curseg-cente))
       or<-order(leveli)
       if (length(or)>1) orseg<-curseg[or,] else orseg<-curseg
       dendat3[curbeg:curend,]<-orseg
       newcolo[curbeg:curend]<-i
       curbeg<-curend+1  #curbeg+length(indit)+1
       i<-i-1
   }

   for (j in 1:d){
       for (i in 1:n){
           x0<-dendat3[i,j]+starts[j]-min(dendat3[,j]) 
               #dendat3[i,j]+origins[j]
           y0<-i
           x1<-origins[j]
           y1<-i
           #segments(x0, y0, x1, y1,col=newcolo[i])
           if (segme)
           segments(x0, y0, x1, y1, col=newcolo[i])
           if (poin)          
           points(x0, y0, col=newcolo[i],pch=pch,cex=cex)

       }
       if (j>1){
            beg<-starts[j]-epsi/2 
            segments(beg,0.5,beg,n+0.5,lty=1,lwd=2)
       }
       segments(origins[j],1,origins[j],n,lty=2)
       if (texto){
       text(origins[j],ystart,as.character(round(means[j],digits=2)),cex=cex)
       text(starts[j]+shift,ystart,as.character(round(minis[j],digits=2)),
       cex=cex)
       } 
   }

}  # else

}






