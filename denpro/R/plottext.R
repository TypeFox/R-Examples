plottext<-function(parents,vecs,lift=0,leimat=NULL,symbo=NULL,cex=NULL){
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc
#
nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)
#
for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}                          
#
#mindiff<-vecs[nodenum,2]-vecs[1,2]
#for (i in 1:(nodenum-1)){
#  diff<-vecs[(i+1),2]-vecs[i,2]
#  if (diff>0) mindiff<-min(diff,mindiff)  
#}
#lift<-mindiff/5
#
moodinum<-length(modloc)
modelocx<-matrix(0,moodinum,1)
modelocy<-matrix(0,moodinum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:moodinum,sep="")
   }
   else{
         if (symbo=="empty") labels<-paste("",1:moodinum,sep="")
         else  labels<-paste(symbo,1:moodinum,sep="")
   }
} 
else{
   labels<-leimat
}
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i] 
    xcor[i]<-xcoor[2*loc-1] 
}
modloc<-omaord2(modloc,xcor)
for (i in 1:moodinum){
    loc<-modloc[i] 
    modelocx[i]<-(xcoor[2*loc-1]+xcoor[2*loc])/2
    modelocy[i]<-ycoor[2*loc-1]+lift
}
text(modelocx,modelocy,labels=labels,cex=cex)
return(list(modelocx=modelocx,labels=labels))
}





