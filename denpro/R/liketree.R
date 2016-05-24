liketree<-function(dendat,pcf,lst)
{

# "lst$infopointer" gives links from nodes to recs
# invert the links
rnum<-length(pcf$value)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

n<-dim(dendat)[1]
d<-dim(dendat)[2]

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# find links from dendat to pcf
# (for simplicity we delete multiple observations in a bin
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
varauslista<-matrix(0,rnum,1)
dendat2<-matrix(0,n,d)
n2<-0
for (i in 1:n){
    j<-1
    notjet<-TRUE
    while ((j<=rnum) && (notjet)){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
              notjet<-FALSE
              if (varauslista[j]==0){ 
                  varauslista[j]<-1
                  n2<-n2+1
                  dendat2[n2,]<-dendat[i,]
                  den2pcf[n2]<-j
                  pcf2den[j]<-n2
              }
         }
         j<-j+1
    }
}
dendat2<-dendat2[1:n2,]

# make tree
parent<-matrix(0,n2,1)
center<-matrix(0,d,n2)
level<-matrix(0,n2,1)
for (i in 1:n2){
   rec<-den2pcf[i]
   node<-nodefinder[rec]
   level[i]<-lst$level[node]

   obs<-0
   curnode<-node
   notfound<-TRUE
   while ((notfound) && (lst$parent[curnode]>0)){
         curnode<-lst$parent[curnode] 
         rec<-lst$infopointer[curnode]
         if (pcf2den[rec]>0){ 
              notfound<-FALSE
              obs<-pcf2den[rec]
         }
   }
   parent[i]<-obs
}
center<-t(dendat2)

return(list(parent=parent,center=center,level=level,
dendat=dendat2,infopoint=seq(1:n2)))
}

