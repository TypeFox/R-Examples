stepwiseR<-function(dendat,leafs,M,pis,mcn,minobs=0,seedi=1,
method="projec",bound=0)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

tr<-myosplitR(dendat,leafs,method,minobs)

suppo<-supp(dendat,blown=TRUE)
step<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(suppo[2*i]-suppo[2*i-1])/n
}

i<-1
while (i<=(M-1)){

   seedi<-seedi+1
   mcdendat<-simutree(tr,mcn,seedi)

   mix<-pis[i]
   trnew<-myosplitpenaR(dendat,leafs,mcdendat,mix,suppo,step,minobs)

   tr<-treeadd(tr,trnew,mix=mix)

   i<-i+1
}

return(tr)

}







