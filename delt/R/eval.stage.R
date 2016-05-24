eval.stage<-function(dendat,leaf,M,pis=NULL,mcn=dim(dendat)[1],
minobs=NULL,seedi=1,
method="projec",bound=0)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

if (is.null(pis)){
   pis<-matrix(0,(M-1),1)
   for (k in 1:(M-1)) pis[k]<-2/(k+2)
}

tr<-eval.greedy(dendat,leaf,method,minobs)

suppo<-supp(dendat,blown=TRUE)
N<-rep(n,d)
step<-stepcalc(suppo,N)

i<-1
while (i<=(M-1)){

   seedi<-seedi+1
   mcdendat<-simutree(tr,mcn,seedi)

   mix<-pis[i]
   trnew<-myosplitpena(dendat,leaf,mcdendat,mix,suppo,step,minobs,method)
   #trnew<-myosplitpenaR(dendat,leaf,mcdendat,mix,suppo,step,minobs,method)

   tr<-treeadd(tr,trnew,mix=mix)

   i<-i+1
}

##################################################
ll<-leaflocs(tr$left,tr$right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)
down<-matrix(0,leafnum,d)
high<-matrix(0,leafnum,d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tr$mean[node]>0){
     efek<-efek+1

     value[efek]<-tr$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tr$low[node,j]
         high[efek,j]<-tr$upp[node,j]
     }
   }
   i<-i+1
}
tr$value<-value[1:efek]
tr$down<-down[1:efek,]
tr$high<-high[1:efek,]
###################################################

tr$N<-N

return(tr)
}








