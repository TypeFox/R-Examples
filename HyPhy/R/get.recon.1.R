get.recon.1 <-
function(phy,phy.sub,reconcile)
{

get.next<-function(prev)
{
if (prev== -(length(phy$tip.label)+1)) return(0)
if (prev >0) return(-phy$edge[prev,1])
which(phy$edge[,2]==-prev)
}

reconcile<-c(reconcile,rep(NA,max(phy.sub$edge)-length(reconcile)))
if (dim(phy.sub$edge)[1]==1)
{
if (is.na(reconcile[1])) reconcile[1]<- -which(phy$tip.label==phy.sub$tip.label[1])
if (is.na(reconcile[2])) reconcile[2]<- 0
return(reconcile)
}
for (i in dim(phy.sub$edge)[1]:0)
{
if (i==0) node<-length(phy.sub$tip.label)+1
else node<-phy.sub$edge[i,2]
if (!is.na(reconcile[node])) next
if (node <= length(phy.sub$tip.label))
reconcile[node]<- -which(phy$tip.label==phy.sub$tip.label[node])
else
{
next.rec<-reconcile[phy.sub$edge[phy.sub$edge[,1]==node,2]]
if (next.rec[1]<0) chain<-get.next(next.rec[1])
else chain<-next.rec[1]
while (chain[1]!=0) chain<-c(get.next(chain[1]),chain)
if (next.rec[2]<0) reconcile[node]<-get.next(next.rec[2])
else reconcile[node]<-next.rec[2]
while (!any(chain==reconcile[node])) 
reconcile[node]<-get.next(reconcile[node])
}
}
reconcile
}
