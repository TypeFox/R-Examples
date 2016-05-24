recon.score <-
function(phy,phy.sub,reconcile=NULL)
{
if(class(phy.sub)=="multiPhylo")
{
out<-NULL
for (i in 1:length(phy.sub)) out<-rbind(out,recon.score(phy,phy.sub[[i]],reconcile))
return(out)
}

if (!is.null(reconcile)) reconcile<- -reconcile

out<-rep(0,2)
reconcile<-get.recon.1(phy,phy.sub,reconcile)
root<-length(phy$tip.label)+1
for(i in dim(phy.sub$edge)[1]:1)
{
ends<-reconcile[phy.sub$edge[i,]]

if (ends[2]>=0) out[1]<-out[1]+1
else if (ends[2]==-root) ends[2]<-0
else if (ends[2]<0) ends[2]<- which(phy$edge[,2]== -ends[2])

while(ends[1]!=ends[2])
{

if (ends[2]==0) stop("I reached the root without finding previous node!") 
else if (ends[2]>0) ends[2]<- -phy$edge[ends[2],1]
else if (ends[2]<0)
{
if (ends[2]==-root) ends[2]<-0
else ends[2]<- which(phy$edge[,2]== -ends[2])
out[2]<-out[2]+1
}
}
}
if (reconcile[length(phy.sub$tip.label)+1]>=0) out[1]<-out[1]+1
out
}
