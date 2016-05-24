kereva<-function(dendat,h,N,kernel="epane",trunc=3,threshold=0.0000001,
hw=NULL,weig=NULL)
{
#weig=rep(1/dim(dendat)[1],dim(dendat)[1]))

#source("~/kerle/profkernCRC.R")
#dyn.load("/home/jsk/kerle/kerCeva")
#dyn.load("/home/jsk/kerle/kerleCversio")
#pk2<-profkernCRC(dendat,h,N,Q)

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3

n<-dim(dendat)[1]
d<-dim(dendat)[2]  #length(N)

if (kernel=="gauss") h<-h*trunc   #trunc<-3

if (is.null(weig)) weig<-rep(1/n,n) 

if (!is.null(hw)){
   weig<-weightsit(n,hw)

   dendatnew<-dendat
   weignew<-weig
   cumul<-0
   for (i in 1:n){
        if (weig[i]>0){
            cumul<-cumul+1
            dendatnew[cumul,]<-dendat[i,]
            weignew[cumul]<-weig[i] 
        }
   }
   dendat<-dendatnew[1:cumul,]
   weig<-weignew[1:cumul]
   n<-cumul
}

inweig<-matrix(0,n+1,1)
inweig[2:(n+1)]<-weig

hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive
{
if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

if (kernel=="radon") kertype<-3
else if (kernel=="epane") kertype<-1 
else kertype<-2  # gaussian

kg<-.C("kergrid",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d),
               as.integer(kertype),
               as.double(trunc),
               as.double(threshold),
               as.double(inweig),
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1),
PACKAGE = "denpro")

#left<-kg$ioleft[2:(kg$numnode+1)]
#right<-kg$ioright[2:(kg$numnode+1)]
#parent<-kg$ioparent[2:(kg$numnode+1)]
#infopointer<-kg$infopointer[2:(kg$numnode+1)]
#iolow<-kg$iolow[2:(kg$numnode+1)]
#ioupp<-kg$ioupp[2:(kg$numnode+1)]

value<-kg$value[2:(kg$numpositive+1)]
#nodefinder<-kg$nodefinder[2:(kg$numpositive+1)]
vecindex<-kg$index[2:(d*kg$numpositive+1)]
index<-matrix(0,kg$numpositive,d)
for (i in 1:kg$numpositive){
  for (j in 1:d){
     index[i,j]<-vecindex[(i-1)*d+j]
  }
}

#return(list(left=left,right=right,parent=parent,infopointer=infopointer,
#low=low,upp=upp,value=value,index=index,nodefinder=nodefinder))

suppo<-matrix(0,2*d,1)
for (i in 1:d){
   suppo[2*i-1]<-min(dendat[,i])-h
   suppo[2*i]<-max(dendat[,i])+h
}

step<-matrix(0,d,1)
for (i in 1:d) step[i]=(suppo[2*i]-suppo[2*i-1])/N[i];

recnum<-dim(index)[1]
low<-matrix(0,recnum,d)
upp<-matrix(0,recnum,d)
for (i in 1:recnum){
     low[i,]<-index[i,]-1
     upp[i,]<-index[i,]
}

return(list(value=value,index=index,
down=low,high=upp,N=N,step=step,support=suppo,n=n))

}
