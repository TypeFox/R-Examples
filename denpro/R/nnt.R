nnt<-function(dendat,f)
{
# dendat is n*d
# f is n-vector of evaluations of the function at dendat

n<-dim(dendat)[1]
d<-dim(dendat)[2]
vol<-pi^(d/2)/gamma(d/2+1)
lkm<-n
parent<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
#center<-matrix(0,d,lkm)
#level<-matrix(0,lkm,1)
levset.radius<-matrix(0,lkm,1)

ford<-order(f) # indeksit pienimmasta suurimpaan
root<-ford[n]
leaf<-root
nindit<-nn.indit(dendat)
neig<-nindit[root,]
notvisited<-setdiff(seq(1,lkm),root)  # a[!a %in% b]  setdiff(a, b) 
radi<-sqrt(sum((dendat[neig[1],]-dendat[root,])^2))/2
volume[root]<-vol*radi^d
levset.radius[root]<-radi

cur<-1
for (i in 1:(n-1)){
    smaller<-ford[n-i]
    nearest<-neig[cur]
    if (smaller==nearest){
        parent[root]<-smaller
        dist.to.parent<-sqrt(sum((dendat[smaller,]-dendat[root,])^2))
        radi<-dist.to.parent+levset.radius[root]
        volume[smaller]<-max(vol*radi^d,volume[smaller])
        levset.radius[smaller]<-max(radi,levset.radius[smaller])
        #  pi*dist(rbind(dendat[root,],dendat[smaller,]))^2
        notvisited<-setdiff(notvisited,root) 
        root<-smaller
        cur<-cur+1
    }
    else{
        parent[root]<-nearest
        #volume[root]<-pi*sum((dendat[root,]-dendat[nearest,])^2)
        notvisited<-setdiff(notvisited,root)
        dist.to.parent<-sqrt(sum((dendat[nearest,]-dendat[root,])^2))
        radi<-dist.to.parent+levset.radius[root]
        volume[nearest]<-max(vol*radi^d,volume[nearest])
        levset.radius[nearest]<-max(radi,levset.radius[nearest])
# dist.nearest.to.root.bound<-2*sqrt(sum((dendat[root,]-dendat[nearest,])^2))
# newvolume<-pi*dist.nearest.to.root.bound^2
# volume[nearest]<-max(volume[nearest],newvolume)   
        root<-smaller
        leaf<-root
        visited<-setdiff(seq(1,lkm),notvisited)
        neig<-setdiff(nindit[root,],visited) #nindit[root,]
        if (volume[root]==0){
           radi<-sqrt(sum((dendat[neig[1],]-dendat[root,])^2))/2
           volume[root]<-vol*radi^d
           levset.radius[root]<-radi
        }
        cur<-1 
    }
}

#mt<-multitree(parent)
#for (i in 1:lkm){
#    if (mt$sibling[i]=0){
#       node<-parent[i]
#       volume[node]

lf<-list(
parent=parent,volume=volume,center=t(dendat),level=f,
#root=root,
#infopointer=infopointer,
#refe=refe,maxdis=maxdis,
dendat=dendat)

return(lf)
}

