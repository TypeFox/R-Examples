requiprobable <-
function(n,labels,label.counts=rep(1,length(labels)))
{

if (all(label.counts<2)) 
{
labels<-labels[label.counts==1]
out<-list()
for (i in 1:n)
{
edge<-matrix(c(0,-1,-1,-1,1,2),3,2)
down<- -2
if (length(labels)>2) for (j in 3:length(labels))
{
addto<-sample(dim(edge)[1],1)
above<-edge[1:addto,]
below<-edge[addto:(dim(edge)[1]),]
above[length(above)]<-below[1]<-down
edge<-rbind(above,c(down,j),below)
down<-down-1
}
edge<-edge[-1,]

up<-length(labels)+1
while (any(edge<0))
{
edge[edge==edge[which(edge<0)[1]]]<-up
up<-up+1
}
phy<-list(edge=edge,tip.label=labels,Nnode=length(labels)-1)
class(phy)="phylo"

   out[[i]]<-phy
   }
   if (n==1) return(out[[1]])
   class(out)<-"multiPhylo"
   return(out)
}


temp<-order(label.counts,decreasing =T)
labels<-labels[temp]
label.counts<-label.counts[temp]

counts<-array(0,dim=label.counts+1)
counts[2]<-1
indices<-list()
for (j in 1:length(label.counts)) indices[[j]]<-1
names(indices)<-1:length(label.counts)

pos<-matrix(1,1,length(label.counts))
pos[1]<-2
while(any(pos!=dim(counts)))
{

pos[1]<-pos[1]+1

if (length(pos)>1)
for(i in 1:(length(label.counts)-1))
if (pos[i]>dim(counts)[i])
{
pos[i]<-1
pos[i+1]<-pos[i+1]+1
}

if (length(pos)>1 && any(pos[-length(pos)]<pos[-1])) 
counts[pos]<-counts[matrix(sort(pos,T),1)]
else
{
for (i in 1:length(label.counts)) indices[[i]]<-1:pos[i]
arr.temp<-extract(counts,indices=indices)

counts[pos]<-sum(arr.temp*arr.temp[length(arr.temp):1])/2

if (all((pos-1)/2==floor((pos-1)/2)))
counts[pos]<-counts[pos]+counts[(pos-1)/2+1]/2


}
}

make.tree<-function(rootn,pos1)
{
if (sum(pos1-1)==1) return(c(NA,which(pos1==2)))

for (i in 1:length(label.counts)) indices[[i]]<-1:pos1[i]
arr.temp<-extract(counts,indices=indices)

probs<-arr.temp*arr.temp[length(arr.temp):1]

if (all((pos1-1)/2==floor((pos1-1)/2)))
probs[(pos1+1)/2]<-probs[(pos1+1)/2]+counts[(pos1+1)/2]

use<-NULL
for (j in 1:length(pos1))
{
choose<-cumsum(apply(probs,j,sum))
use<-c(use,which(choose>runif(1,0,choose[length(choose)]))[1])
indices[[j]]<-use[length(use)]
probs<-extract(probs,indices=indices)
indices[[j]]<-1
}

x1<-make.tree(rootn+1,use)
x2<-make.tree(rootn+sum(use-1),pos1+1-use)

x1[1]<-x2[1]<-rootn

rbind(c(NA,rootn),x1,x2)
}

   out<-list()
   for (k in 1:n)
   {
edge<-make.tree(sum(label.counts)+1,matrix(label.counts,1)+1)[-1,]
tip.label<-NULL

if (length(label.counts)>1)
for (i in length(label.counts):2)
edge[which(edge[,2]==i),2]<-1:label.counts[i]+sum(label.counts[1:(i-1)])
edge[which(edge[,2]==1),2]<-1:label.counts[1]

tip.label<-NULL

for (i in 1:length(label.counts)) tip.label<-c(tip.label,rep(labels[i],label.counts[i]))

phy<-list(edge=edge,tip.label=tip.label,Nnode=sum(label.counts)-1)
class(phy)="phylo"

   out[[k]]<-phy
   }
   if (n==1) return(out[[1]])
   class(out)<-"multiPhylo"
   out
}
