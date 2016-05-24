#LOOP ANALYSIS
#version 1.1
#2012/8/10
#YOUHUA CHEN
###############

library(grid)
library(MASS)

#forward algorithm to search loops
#always search the first descendant firstly, then other descandants
loop.forward<-function(gemat)
{
nodenum<-dim(gemat)[1]
series<-1:nodenum
for(i in series)
{
if(gemat[i,i]>0)
{
thisloop<-c(i,i)
weight<-gemat[i,i]
thismat<-matrix(0,ncol=nodenum,nrow=nodenum)
thismat[i,i]=weight
return(list(loop=thisloop,weights=weight,loopmat=thismat,newmat=gemat-thismat))
break
}else
{
unchecked<-series
thisnode<-i
father<-thisnode
out.index<-which(gemat[thisnode,]>0)
out.index<-out.index[which(out.index!=thisnode)] #delete itself
in.index<-which(gemat[,thisnode]>0)
in.index<-in.index[which(in.index!=thisnode)] #delete itself
thisloop<-vector()
weight<-vector()
flist<-list() #desendant and anscestor list
if(length(out.index)>0 & length(in.index)>0)
{
thisloop<-c(thisloop,thisnode)
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
#dynamic search
thisnode<-flist[[length(flist)]]$desc[1]
father=thisnode
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1] #delete the first desendant iteratively
#begin to search other nodes
flag=1
while(gemat[thisnode,i]==0 & flag==1)#not the last node
{
out.index<-which(gemat[thisnode,]>0)
out.index<-out.index[which(out.index!=thisnode)] #delete itself
if(length(out.index)>0 & thisnode %in% unchecked==TRUE)
{
#dynamic search
thisloop<-c(thisloop,thisnode)
unchecked<-unchecked[-thisnode]
weight<-c(weight,gemat[thisloop[length(thisloop)-1],thisloop[length(thisloop)]])
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
thisnode<-flist[[length(flist)]]$desc[1]
father=thisnode
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1] #delete the first desendant iteratively
flag=1
#
}else
{
#dynamic search
#backward search
backone<-length(flist)
while(length(flist[[backone]]$desc)==0)
{
flist[[backone]]<-NULL
thisloop<-thisloop[-length(thisloop)]
weight<-weight[-length(weight)]
backone<-backone-1
if(backone==0)
{
flag=0
break
}
}#while backone
if(flag==1)
{
thisnode<-flist[[backone]]$desc[1]
father=thisnode
flist[[backone]]$desc<-flist[[backone]]$desc[-1] #delete the first desendant iteratively
}
#
}#if else
#
if(length(flist)>0)
{
flag=1
}else
{
flag=0
break
}
}#while
if(flag==1)
{
thisloop<-c(thisloop,thisnode,i)
weight<-c(weight,gemat[thisloop[length(thisloop)-2],thisnode])
weight<-c(weight,gemat[thisnode,i])
thismat<-matrix(0,ncol=nodenum,nrow=nodenum)
for(i in 1:(length(thisloop)-1))
{
thismat[thisloop[i],thisloop[i+1]]=min(weight)
}
return(list(loop=thisloop,weights=weight,loopmat=thismat,newmat=gemat-thismat))
break
}#no return
}#if
}#not self-looped
}#i
}#end function

#random algorithm to search loops
#search one descendant randomly each time
loop.random<-function(gemat)
{
nodenum<-dim(gemat)[1]
series<-1:nodenum
for(i in series)
{
if(gemat[i,i]>0)
{
thisloop<-c(i,i)
weight<-gemat[i,i]
thismat<-matrix(0,ncol=nodenum,nrow=nodenum)
thismat[i,i]=weight
return(list(loop=thisloop,weights=weight,loopmat=thismat,newmat=gemat-thismat))
break
}else
{
unchecked<-series
thisnode<-i
father<-thisnode
out.index<-which(gemat[thisnode,]>0)
out.index<-out.index[which(out.index!=thisnode)] #delete itself
in.index<-which(gemat[,thisnode]>0)
in.index<-in.index[which(in.index!=thisnode)] #delete itself
thisloop<-vector()
weight<-vector()
flist<-list() #desendant and anscestor list
if(length(out.index)>0 & length(in.index)>0)
{
thisloop<-c(thisloop,thisnode)
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
#dynamic search
rid<-sample(1:length(out.index),1,replace=FALSE)
thisnode<-flist[[length(flist)]]$desc[rid]
father=thisnode
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-rid] #delete the first desendant iteratively
#begin to search other nodes
flag=1
while(gemat[thisnode,i]==0 & flag==1)#not the last node
{
out.index<-which(gemat[thisnode,]>0)
out.index<-out.index[which(out.index!=thisnode)] #delete itself
if(length(out.index)>0 & thisnode %in% unchecked==TRUE)
{
#dynamic search
thisloop<-c(thisloop,thisnode)
unchecked<-unchecked[-thisnode]
weight<-c(weight,gemat[thisloop[length(thisloop)-1],thisloop[length(thisloop)]])
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
rid<-sample(1:length(out.index),1,replace=FALSE)
thisnode<-flist[[length(flist)]]$desc[rid]
father=thisnode
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-rid] #delete the first desendant iteratively
flag=1
#
}else
{
#dynamic search
#backward search
backone<-length(flist)
while(length(flist[[backone]]$desc)==0)
{
flist[[backone]]<-NULL
thisloop<-thisloop[-length(thisloop)]
weight<-weight[-length(weight)]
backone<-backone-1
if(backone==0)
{
flag=0
break
}
}#while backone
if(flag==1)
{
rid<-sample(1:length(out.index),1,replace=FALSE)
thisnode<-flist[[backone]]$desc[rid]
father=thisnode
flist[[backone]]$desc<-flist[[backone]]$desc[-rid] #delete the first desendant iteratively
}
#
}#if else
#
if(length(flist)>0)
{
flag=1
}else
{
flag=0
break
}
}#while
if(flag==1)
{
thisloop<-c(thisloop,thisnode,i)
weight<-c(weight,gemat[thisloop[length(thisloop)-2],thisnode])
weight<-c(weight,gemat[thisnode,i])
thismat<-matrix(0,ncol=nodenum,nrow=nodenum)
for(i in 1:(length(thisloop)-1))
{
thismat[thisloop[i],thisloop[i+1]]=min(weight)
}
return(list(loop=thisloop,weights=weight,loopmat=thismat,newmat=gemat-thismat))
break
}#no return
}#if
}#not self-looped
}#i
}#end function


#decomposition based on random search method
decomp<-function(gemat)
{
results<-list()
remainder<-list()
loops<-list()
i=0
results<-loop.random(gemat)
while(is.list(results)==TRUE)
{
i=i+1
gemat<-results$newmat
loops[[i]]<-results$loopmat
remainder<-results$newmat
results<-loop.random(gemat)
}
loops[[length(loops)+1]]<-remainder
return(loops)
}

#rank nodes based on links number
rank.nodes<-function(gemat,type="both")
{
num<-dim(gemat)[1]
nodes<-vector()
for(i in 1:num)
{
if(gemat[i,i]>0)
{
dd=1
}else
{
dd=0
}
inum<-length(which(gemat[,i]>0))
onum<-length(which(gemat[i,]>0))
if(type=="both")
{
nodes[i]<-inum+onum-dd*2
}
if(type=="in")
{
nodes[i]<-inum-dd
}
if(type=="out")
{
nodes[i]<-onum-dd
}
}#i
names(nodes)<-seq(1,num,1)
return(sort(nodes,decreasing=TRUE))
}

#plot directed graph in r from edge graph matrix
gplot<-function(edgemat,arrow=TRUE,lty=1,col=8,weighted=TRUE)
{
grid.newpage()
pi=3.141593
mat<-edgemat
mat[,1]<-edgemat[,2]
mat[,2]<-edgemat[,1]
num<-unique(as.vector(mat[,1:2]))
edges<-dim(mat)[1]
if(weighted==TRUE)
{
w<-mat[,3]/max(mat[,3])*5 #weights in each edge
}else
{
w<-rep(1,1,length(num))
}
#
if(length(num)>50)
{
r<-.4
arrowsize<-.001
lim<-.05
}else
{
r<-.3
arrowsize<-.02
lim<-.2
}
#check self-loop
self=FALSE
for(i in 1:edges)
{
if(mat[i,1]==mat[i,2])
{
self=TRUE
break
}
}

#create the coordinates for each species
xy<-matrix(0,ncol=2,nrow=max(num))
pxy<-xy
name<-xy
t1<-xy
t2<-xy
o<-xy
for(i in 1:length(num))
{
ang<-2*pi/length(num)*(i-1)
xy[num[i],1]<-.5+r*sin(ang)
xy[num[i],2]<-.5+r*cos(ang)
if(self==FALSE)
{
name[num[i],1]<-.5+(r+.05)*sin(ang)
name[num[i],2]<-.5+(r+.05)*cos(ang)
}else
{
name[num[i],1]<-.5+(r+.1)*sin(ang)
name[num[i],2]<-.5+(r+.1)*cos(ang)
t1[num[i],1]<-.5+(r+.014)*sin(ang)
t1[num[i],2]<-.5+(r+.014)*cos(ang)
t2[num[i],1]<-.5+(r+.08)*sin(ang)
t2[num[i],2]<-.5+(r+.08)*cos(ang)
o[num[i],1]<-.5+(r+.047)*sin(ang) #radius=.01
o[num[i],2]<-.5+(r+.047)*cos(ang)
radius=.047
}
pxy[num[i],1]<-.5+(r+.01)*sin(ang)
pxy[num[i],2]<-.5+(r+.01)*cos(ang)
grid.points(unit(pxy[num[i],1],"npc"),unit(pxy[num[i],2],"npc"),pch=20,gp=gpar(col=4))
if(name[num[i],1]<=.5)
{
grid.text(as.character(num[i]),unit(name[num[i],1],"npc"),unit(name[num[i],2],"npc"),just="left")
}else
{
grid.text(as.character(num[i]),unit(name[num[i],1],"npc"),unit(name[num[i],2],"npc"),just="right")
}
}
#draw lines or curves from graph matrix
for(i in 1:edges)
{
if(mat[i,1]==mat[i,2])
{
if(t1[mat[i,1],1]==t2[mat[i,1],1])
{
lx<-o[mat[i,1],1]-radius
ly<-o[mat[i,1],2]
rx<-o[mat[i,1],1]+radius
ry<-o[mat[i,1],2]
}else if(t1[mat[i,1],2]==t2[mat[i,1],2])
{
lx<-o[mat[i,1],1]
ly<-o[mat[i,1],2]+radius
rx<-o[mat[i,1],1]
ry<-o[mat[i,1],2]-radius
}else
{
k=(t1[mat[i,1],2]-t2[mat[i,1],2])/(t1[mat[i,1],1]-t2[mat[i,1],1])
lx<-o[mat[i,1],1]-radius*k/(sqrt(k^2+1))
ly<-o[mat[i,1],2]+radius/(sqrt(k^2+1))
rx<-o[mat[i,1],1]+radius*k/(sqrt(k^2+1))
ry<-o[mat[i,1],2]-radius/(sqrt(k^2+1))
}
xx<-c(t1[mat[i,1],1],lx,t2[mat[i,2],1],rx,t1[mat[i,1],1])
yy<-c(t1[mat[i,1],2],ly,t2[mat[i,2],2],ry,t1[mat[i,1],2])
if(arrow==TRUE)
{
grid.xspline(xx,yy,shape=-1,arrow=arrow(angle=30,length=unit(arrowsize,"npc")),gp=gpar(lty=lty,col=col,lwd=w[i]))
}else
{
grid.xspline(xx,yy,shape=-1,gp=gpar(lty=lty,col=col,lwd=w[i]))
}
}else
{
if(length(which(mat[i:edges,1]==mat[i,2] & mat[i:edges,2]==mat[i,1]))==0)
{
if(arrow==TRUE)
{
grid.segments(xy[mat[i,1],1],xy[mat[i,1],2],xy[mat[i,2],1],xy[mat[i,2],2],arrow=arrow(angle=15,length=unit(arrowsize,"npc")),gp=gpar(lty=lty,col=col,lwd=w[i]))
}else
{
grid.segments(xy[mat[i,1],1],xy[mat[i,1],2],xy[mat[i,2],1],xy[mat[i,2],2],gp=gpar(lty=lty,col=col,lwd=w[i]))
}
}else
{
if(arrow==TRUE)
{
grid.xspline(x=c(xy[mat[i,1],1],(xy[mat[i,1],1]+xy[mat[i,2],1])/2,xy[mat[i,2],1]),y=c(xy[mat[i,1],2],(xy[mat[i,1],2]+xy[mat[i,2],2]+runif(1,lim/10,lim))/2,xy[mat[i,2],2]),shape=1,arrow=arrow(angle=15,length=unit(arrowsize,"npc")),gp=gpar(lty=lty,col=col,lwd=w[i]))
}else
{
grid.xspline(x=c(xy[mat[i,1],1],(xy[mat[i,1],1]+xy[mat[i,2],1])/2,xy[mat[i,2],1]),y=c(xy[mat[i,1],2],(xy[mat[i,1],2]+xy[mat[i,2],2]+runif(1,lim/10,lim))/2,xy[mat[i,2],2]),shape=1,gp=gpar(lty=lty,col=col,lwd=w[i]))
}
}
}
}#i
}#end function

#plot directed graphs directly from square graph matrix
gplot1<-function(gemat,arrow=TRUE,lty=1,col=8,weighted=TRUE)
{
mat<-convertion(gemat=gemat)
gplot(edgemat=mat,arrow=arrow,lty=lty,col=col,weighted=weighted)
}

#convert graph matrix form into edge form
convertion<-function(gemat)
{
mat<-gemat
num<-length(which(as.vector(mat)>0))
edges<-matrix(0,nrow=num,ncol=3)
t=0
for(i in 1:dim(mat)[1])
{
for(j in 1:dim(mat)[2])
{
if(mat[i,j]>0)
{
t=t+1
edges[t,1]<-i
edges[t,2]<-j
edges[t,3]<-mat[i,j]
}
}
}
return(edges)
}

#list all pathways in the network: including circle and self-looped pathways
pathways<-function(gemat,bsp)
{
nodenum<-dim(gemat)[1]
series<-1:nodenum
ps<-list()
pw<-list()
for(i in 1:length(bsp))
{
#
thisloop<-vector()
weight<-vector()
flist<-list() #desendant and anscestor list
flag=1
unchecked<-series
thisnode<-bsp[i]
father<-thisnode
out.index<-which(gemat[thisnode,]>0)
#out.index<-out.index[which(out.index!=thisnode)] #delete itself
if(length(out.index)>0)
{
thisloop<-c(thisloop,thisnode)
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
#dynamic search
thisnode<-flist[[length(flist)]]$desc[1]
#
# while(length(which(thisloop==thisnode))>0) #self-loop cases
# {
# ps[[length(ps)+1]]<-c(thisloop,thisnode)
# pw[[length(pw)+1]]<-c(weight,gemat[thisnode,thisnode])
# flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1] #delete the first desendant iteratively
# thisnode<-flist[[length(flist)]]$desc[1]
# }
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1] #delete the first desendant iteratively
father=thisnode
flag=1
#begin to search other nodes
#loop is not possible to larger than nodenum!
while(flag==1)#not the last node
{
out.index<-which(gemat[thisnode,]>0)
#out.index<-out.index[which(out.index!=thisnode)] #delete itself
if(length(out.index)>0)
{
#dynamic search
thisloop<-c(thisloop,thisnode)
weight<-c(weight,gemat[thisloop[length(thisloop)-1],thisloop[length(thisloop)]])
flist[[length(flist)+1]]<-list(ansc=father,desc=out.index)
thisnode<-flist[[length(flist)]]$desc[1]
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1] #delete the first desendant iteratively
father=thisnode
flag=1
#
}else
{
#dynamic search
#backward search
backone<-length(flist)
onetime<-backone
#
###
#for the case when $desc!=0, but out.index==0
ps[[length(ps)+1]]<-c(thisloop,thisnode)
pw[[length(pw)+1]]<-c(weight,gemat[thisloop[length(thisloop)],thisnode])
###
#
while(length(flist[[backone]]$desc)==0)
{
# if(onetime>0 & is.na(thisnode)==FALSE)
# {
# ps[[length(ps)+1]]<-c(thisloop,thisnode)
# pw[[length(pw)+1]]<-c(weight,gemat[thisloop[length(thisloop)],thisnode])
# }
#
unchecked<-unchecked[-flist[[backone]]$ansc]
flist[[backone]]<-NULL
onetime<-0
thisloop<-thisloop[-length(thisloop)]
weight<-weight[-length(weight)]
backone<-backone-1
if(backone==0)
{
flag=0
break
}
}#while backone
#
if(flag==1)
{
thisnode<-flist[[backone]]$desc[1]
father=thisnode
flist[[backone]]$desc<-flist[[backone]]$desc[-1] #delete the first desendant iteratively
}
#
}#if else
#
if(length(flist)>0)
{
flag=1
}else
{
flag=0
}
########################
#handling special cases
#avoid chain self-loop
#back check of self-loop
if(length(which(thisloop[-length(thisloop)]==thisloop[length(thisloop)]))>0)
{
ps[[length(ps)+1]]<-thisloop
pw[[length(pw)+1]]<-weight
thisloop<-thisloop[-length(thisloop)]
weight<-weight[-length(weight)]
flist[[length(flist)]]<-NULL
thisnode<-flist[[length(flist)]]$desc[1]
flist[[length(flist)]]$desc<-flist[[length(flist)]]$desc[-1]
father=thisnode
flag=1
}
}#while
#
}#if
#
}#i
return(list(pathways=ps,pathwayweights=pw))
}

#unique pathways in the network for a specific species
uniquepaths<-function(gemat,sp)
{
paths<-pathways(gemat=gemat,bsp=sp)
num<-length(paths$pathways)
up<-list()
upw<-list()
#
up[[1]]<-paths$pathways[[1]]
upw[[1]]<-paths$pathwayweights[[1]]
if(num>1)
{
for(i in 2:num)
{
can<-paths$pathways[[i]]
flag=1
for(j in 1:(i-1))
{
if(length(can)==length(paths$pathways[[j]]))
{
if(length(can[can %in% paths$pathways[[j]]])==length(can))
{
flag=0
break
}
}
}#j
#add into the unique list
if(flag==1)
{
up[[length(up)+1]]<-can
upw[[length(upw)+1]]<-paths$pathwayweights[[i]]
}
}#i
}#if
return(list(paths=up,weights=upw))
}



#find the longest pathways based on the species given
longest.chain<-function(gemat,sp)
{
paths<-uniquepaths(gemat=gemat,sp=sp)
pnum<-length(paths$paths)
svec<-vector(length=pnum)
for(i in 1:pnum)
{
svec[i]<-length(paths$paths[[i]])
}
maxi<-max(svec)
ind<-which(svec==maxi)
num<-length(ind)
candidates<-list()
cweights<-list()
for(i in 1:num)
{
candidates[[i]]<-paths$paths[[ind[i]]]
cweights[[i]]<-paths$weights[[ind[i]]]
}
return(list(paths=candidates,weights=cweights))
}


#find the shortest pathways based on the species given,not including self-loop cases
shortest.chain<-function(gemat,sp)
{
paths<-uniquepaths(gemat=gemat,sp=sp)
pnum<-length(paths$paths)
svec<-vector(length=pnum)
for(i in 1:pnum)
{
this<-length(paths$paths[[i]])
if(this==2)
{
svec[i]=999
}else
{
svec[i]<-this
}
}
maxi<-min(svec)
ind<-which(svec==maxi)
num<-length(ind)
candidates<-list()
cweights<-list()
for(i in 1:num)
{
candidates[[i]]<-paths$paths[[ind[i]]]
cweights[[i]]<-paths$weights[[ind[i]]]
}
return(list(paths=candidates,weights=cweights))
}


#find the pathways with largest weights
largest.weight<-function(gemat,sp)
{
paths<-uniquepaths(gemat=gemat,sp=sp)
pnum<-length(paths$weights)
svec<-vector(length=pnum)
for(i in 1:pnum)
{
svec[i]<-sum(paths$weights[[i]])
}
maxi<-max(svec)
ind<-which(svec==maxi)
num<-length(ind)
candidates<-list()
cweights<-list()
for(i in 1:num)
{
candidates[[i]]<-paths$paths[[ind[i]]]
cweights[[i]]<-paths$weights[[ind[i]]]
}
return(list(paths=candidates,weights=cweights))
}


#find the pathways with smallest weights
smallest.weight<-function(gemat,sp)
{
paths<-uniquepaths(gemat=gemat,sp=sp)
pnum<-length(paths$weights)
svec<-vector(length=pnum)
for(i in 1:pnum)
{
svec[i]<-sum(paths$weights[[i]])
}
maxi<-min(svec)
ind<-which(svec==maxi)
num<-length(ind)
candidates<-list()
cweights<-list()
for(i in 1:num)
{
candidates[[i]]<-paths$paths[[ind[i]]]
cweights[[i]]<-paths$weights[[ind[i]]]
}
return(list(paths=candidates,weights=cweights))
}


#find the largest weighted ones among the longest pathways based on the species given
lclw<-function(gemat,sp)
{
paths<-uniquepaths(gemat=gemat,sp=sp)
pnum<-length(paths$paths)
svec<-vector(length=pnum)
for(i in 1:pnum)
{
svec[i]<-length(paths$paths[[i]])
}
maxi<-max(svec)
ind<-which(svec==maxi)
num<-length(ind)
candidates<-list()
cweights<-list()
wsum<-vector()
for(i in 1:num)
{
candidates[[i]]<-paths$paths[[ind[i]]]
cweights[[i]]<-paths$weights[[ind[i]]]
wsum[i]<-sum(cweights[[i]])
}
ind2<-which(wsum==max(wsum))
c2<-list()
w2<-list()
for(i in 1:length(ind2))
{
c2[[i]]<-candidates[[ind2[i]]]
w2[[i]]<-cweights[[ind2[i]]]
}
return(list(chains=c2,weights=w2))
}


#Primm algorithm to find minimum spanning tree
mst.primm<-function(gemat)
{
used<-vector()
mstmat<-vector()
nnum<-dim(gemat)[1]
notused<-seq(1,nnum,1)
this<-1
used<-c(used,this)
notused<-notused[-this]
while(length(used)<nnum)
{
conmat<-as.matrix(gemat[used,notused])
if(length(used)==1)
{
conmat<-t(conmat)
}
minimum<-min(as.vector(conmat[which(conmat>0)]))
for(i in 1:length(used))
{
if(length(which(conmat[i,]==minimum))>0)
{
x<-used[i]
y<-notused[which(conmat[i,]==minimum)]
break
}
}#i
mstmat<-rbind(mstmat,c(x,y))
used<-c(used,y)
notused<-notused[which(notused!=y)]
}#while
return(mstmat)
}#end function


#node similarity based on pointed nodes' similarity
node.similarity<-function(gemat,type="both",metric="jaccard")
{
mat<-convertion(gemat=gemat)
num<-dim(gemat)[1]
sim<-matrix(0,num,num)
for(i in 1:(num-1))
{
for(j in (i+1):num)
{
if(type=="both")
{
ind1<-which(mat[,1]==i)
ind2<-which(mat[,2]==i)
jnd1<-which(mat[,1]==j)
jnd2<-which(mat[,2]==j)
ii1<-mat[ind1,2]
jj1<-mat[jnd1,2]
ii2<-mat[ind2,1]
jj2<-mat[jnd2,1]
ii<-unique(c(ii1,ii2))
jj<-unique(c(jj1,jj2))
shared<-length(ii[ii %in% jj])
ionly<-length(ii[!ii %in% jj])
jonly<-length(jj[!jj %in% ii])
if(metric=="jaccard")
{
sim[i,j]=(shared)/(ionly+jonly+shared)
}
if(metric=="sorensen")
{
sim[i,j]=(2*shared)/(ionly+jonly+2*shared)
}
}#type=both
if(type=="in")
{
ind<-which(mat[,2]==i)
jnd<-which(mat[,2]==j)
ii<-mat[ind,1]
jj<-mat[jnd,1]
shared<-length(ii[ii %in% jj])
ionly<-length(ii[!ii %in% jj])
jonly<-length(jj[!jj %in% ii])
if(metric=="jaccard")
{
sim[i,j]=(shared)/(ionly+jonly+shared)
}
if(metric=="sorensen")
{
sim[i,j]=(2*shared)/(ionly+jonly+2*shared)
}
}#type=in
if(type=="out")
{
ind<-which(mat[,1]==i)
jnd<-which(mat[,1]==j)
ii<-mat[ind,2]
jj<-mat[jnd,2]
shared<-length(ii[ii %in% jj])
ionly<-length(ii[!ii %in% jj])
jonly<-length(jj[!jj %in% ii])
if(metric=="jaccard")
{
sim[i,j]=(shared)/(ionly+jonly+shared)
}
if(metric=="sorensen")
{
sim[i,j]=(2*shared)/(ionly+jonly+2*shared)
}
}#type=out
}#j
}#i
sim=(sim+t(sim))/2
sim=sim+diag(num)
return(sim)
}


#make ordination of the nodes based on node similarity using NMDS method
#return coordinates of the nodes and their labels
nmds.ordination<-function(gemat,type="both",metric="jaccard")
{
sim<-node.similarity(gemat=gemat,type=type,metric=metric)
num<-dim(gemat)[1]
sim<-matrix(1,num,num)-sim
res<-isoMDS(sim)
xy<-res$points
#
if(length(colnames(gemat))>0)
{
names<-colnames(gemat)
}else
{
names<-seq(1,num,1)
}
return(list(coord=xy,names=names))
}

#plot network using Non-dimensional scaling technique to allow similar nodes
#gather together, while dissimilar nodes separated
#require(MASS)
fplot<-function(gemat,type="both",metric="jaccard",addlabels=FALSE,scaled=TRUE,weighted=TRUE,pch=20,bg=1,pcex=3,pcol=4,lty=1,lcol=8,tfont=12,tcol=1)
{
dat<-nmds.ordination(gemat=gemat,type=type,metric=metric)
num<-dim(gemat)[1]
dat$coord<-dat$coord+matrix(rnorm(num*2,mean=0,sd=.1),ncol=2,nrow=num)
edges<-convertion(gemat=gemat)
plot(dat$coord,,type="n",axes=FALSE,xlab="",ylab="",pch=pch,cex=pcex)
#
maxi<-max(edges[,3])
for(i in 1:dim(edges)[1])
{
x1<-dat$coord[edges[i,1],1]
y1<-dat$coord[edges[i,1],2]
x2<-dat$coord[edges[i,2],1]
y2<-dat$coord[edges[i,2],2]
if(weighted==TRUE)
{
wei<-edges[i,3]/maxi*5
}else
{
wei<-rep(1,1,length(edges[,3]))
}
if(scaled==TRUE)
{
segments(x1,y1,x2,y2,col=lcol,lty=lty,lwd=wei)
}else
{
segments(x1,y1,x2,y2,col=lcol,lty=lty)
}
}#i
points(dat$coord,pch=pch,cex=pcex,col=pcol,bg=bg)
if(addlabels==TRUE)
{
text(dat$coord,labels=dat$names,font=tfont,col=tcol)
}
}


#group network plots
#allow you to plot different networks together
groupplot<-function(gemat,groups,type="both",metric="jaccard",addlabels=FALSE,scaled=TRUE,pch=20,bg=1,pcex=3,pcol=4,lty=1,lcol=8,tfont=12,tcol=1)
{
data<-nmds.ordination(gemat=gemat,type=type,metric=metric)
num<-dim(gemat)[1]
data$coord<-data$coord+matrix(rnorm(num*2,mean=0,sd=.1),ncol=2,nrow=num)
alledges<-convertion(gemat=gemat)
plot(data$coord,type="n",axes=FALSE,xlab="",ylab="",pch=pch,cex=pcex)
#
coords<-data$coord
maxi<-max(alledges[,3])
#
#check lazy parameters!!!
if(length(pch)!=length(groups))
{
this<-pch[1]
pch<-rep(this,1,length(groups))
}
if(length(bg)!=length(groups))
{
this<-bg[1]
bg<-rep(this,1,length(groups))
}
if(length(pcex)!=length(groups))
{
this<-pcex[1]
pcex<-rep(this,1,length(groups))
}
if(length(pcol)!=length(groups))
{
this<-pcol[1]
pcol<-rep(this,1,length(groups))
}
if(length(lty)!=length(groups))
{
this<-lty[1]
lty<-rep(this,1,length(groups))
}
if(length(lcol)!=length(groups))
{
this<-lcol[1]
lcol<-rep(this,1,length(groups))
}
if(length(tfont)!=length(groups))
{
this<-tfont[1]
tfont<-rep(this,1,length(groups))
}
if(length(tcol)!=length(groups))
{
this<-tcol[1]
tcol<-rep(this,1,length(groups))
}
#
for(j in 1:length(groups))
{
sp<-groups[[j]]
edges<-vector()
dat<-coords[sp,]
names<-data$names[sp]
for(i in 1:dim(alledges)[1])
{
this<-alledges[i,1:2]
if(length(this[this %in% sp])==length(this))
{
edges<-rbind(edges,alledges[i,])
}
}#check out all needed links
for(i in 1:dim(edges)[1])
{
x1<-dat[which(names==edges[i,1]),1]
y1<-dat[which(names==edges[i,1]),2]
x2<-dat[which(names==edges[i,2]),1]
y2<-dat[which(names==edges[i,2]),2]
wei<-edges[i,3]/maxi*5
if(scaled==TRUE)
{
segments(x1,y1,x2,y2,col=lcol[j],lty=lty[j],lwd=wei)
}else
{
segments(x1,y1,x2,y2,col=lcol[j],lty=lty[j])
}
}#i
#add points and texts
points(dat,pch=pch[j],cex=pcex[j],col=pcol[j],bg=bg[j])
if(addlabels==TRUE)
{
text(dat,labels=names,font=tfont[j],col=tcol[j])
}
}#j
}


#food web plot
#a special plot for food webs following vertical cascades!
fplot.foodweb<-function(gemat,ranks,addlabels=FALSE,scaled=TRUE,weighted=TRUE,pch=20,bg=1,pcex=3,pcol=4,lty=1,lcol=8,tfont=12,tcol=1)
{
dat<-list()
num<-dim(gemat)[1]
dat$coord<-matrix(0,nrow=num,ncol=2)
if(length(colnames(gemat))>0)
{
dat$names<-colnames(gemat)
}else
{
dat$names<-seq(1,num,1)
}
lnum<-length(ranks)
for(i in 1:lnum)
{
sp<-ranks[[i]]
y<-i/lnum
for(j in 1:length(sp))
{
dat$coord[sp[j],1]<-j/length(sp)
dat$coord[sp[j],2]<-y
}
if(length(sp)>1)
{
for(j in 1:length(sp))
{
if(length(which(gemat[sp[j],sp[-j]]>0))>0)
{
ind<-which(gemat[sp[j],sp[-j]]>0)
this<-sp[-j]
this<-this[ind]
ymin<-min(dat$coord[this,2])
dat$coord[sp[j],2]<-ymin-y*.05 #lower than other species in the same rank
}
}#j
}
midpoint<-(max(dat$coord[sp,1])-min(dat$coord[sp,1]))/2
#move towards the centre
distance<-.5-midpoint
dat$coord[sp,1]<-dat$coord[sp,1]+distance
}#i
#
dat$coord[which(dat$coord[,1]>1),1]<-dat$coord[which(dat$coord[,1]>1),1]-1
#
edges<-convertion(gemat=gemat)
plot(dat$coord,,type="n",axes=FALSE,xlab="",ylab="",pch=pch,cex=pcex)
#
maxi<-max(edges[,3])
for(i in 1:dim(edges)[1])
{
x1<-dat$coord[edges[i,1],1]
y1<-dat$coord[edges[i,1],2]
x2<-dat$coord[edges[i,2],1]
y2<-dat$coord[edges[i,2],2]
if(weighted==TRUE)
{
wei<-edges[i,3]/maxi*5
}else
{
wei=1
}
if(scaled==TRUE)
{
segments(x1,y1,x2,y2,col=lcol,lty=lty,lwd=wei)
}else
{
segments(x1,y1,x2,y2,col=lcol,lty=lty)
}
}#i
points(dat$coord,pch=pch,cex=pcex,col=pcol,bg=bg)
if(addlabels==TRUE)
{
text(dat$coord,labels=dat$names,font=tfont,col=tcol)
}
}


#make food trophic ranks for all the species in the matrix for fplot.foodweb function
#gemat is a square matrix
find.ranks<-function(gemat,converted=TRUE)
{
spnum<-dim(gemat)[1]
ranks<-rep(0,1,spnum)
ranks[1]<-1
for(i in 2:spnum)
{
ind<-which(gemat[i,1:(i-1)]>0)
if(length(ind)==0)
{
ranks[i]<-1
}else
{
ranks[i]<-max(ranks[ind])+1
}
}#i
#adjust top species' ranks
ranks[length(ranks)]=max(ranks)
top<-max(ranks)
for(i in 1:(spnum-1))
{
if(length(which(gemat[(i+1):spnum,i]>0))==0)
{
ranks[i]=top
}
}
if(converted==TRUE)
{
maxi<-max(ranks)
this<-ranks
ranks<-list()
for(i in 1:maxi)
{
ranks[[i]]<-which(this==i)
}
}#converted
return(ranks)
}


#food web group plot
#a special group plot for food webs following vertical cascades!
groupplot.foodweb<-function(gemat,ranks,groups,addlabels=FALSE,scaled=TRUE,pch=20,bg=1,pcex=3,pcol=4,lty=1,lcol=8,tfont=12,tcol=1)
{
dat<-list()
num<-dim(gemat)[1]
dat$coord<-matrix(0,nrow=num,ncol=2)
if(length(colnames(gemat))>0)
{
dat$names<-colnames(gemat)
}else
{
dat$names<-seq(1,num,1)
}
lnum<-length(ranks)
for(i in 1:lnum)
{
sp<-ranks[[i]]
y<-i/lnum
for(j in 1:length(sp))
{
dat$coord[sp[j],1]<-j/length(sp)
dat$coord[sp[j],2]<-y
}
for(j in 1:length(sp))
{
if(length(which(gemat[sp[j],sp[-j]]>0))>0)
{
ind<-which(gemat[sp[j],sp[-j]]>0)
this<-sp[-j]
this<-this[ind]
ymin<-min(dat$coord[this,2])
dat$coord[sp[j],2]<-ymin-y*.05 #lower than other species in the same rank
}
}#j
midpoint<-(max(dat$coord[sp,1])-min(dat$coord[sp,1]))/2
#move towards the centre
distance<-.5-midpoint
dat$coord[sp,1]<-dat$coord[sp,1]+distance
}#i
#
dat$coord[which(dat$coord[,1]>1),1]<-dat$coord[which(dat$coord[,1]>1),1]-1
#
#plot now
data<-dat
num<-dim(gemat)[1]
data$coord<-data$coord+matrix(rnorm(num*2,mean=0,sd=.1),ncol=2,nrow=num)
alledges<-convertion(gemat=gemat)
plot(data$coord,type="n",axes=FALSE,xlab="",ylab="",pch=pch,cex=pcex)
#
coords<-data$coord
maxi<-max(alledges[,3])
#
#check lazy parameters!!!
if(length(pch)!=length(groups))
{
this<-pch[1]
pch<-rep(this,1,length(groups))
}
if(length(bg)!=length(groups))
{
this<-bg[1]
bg<-rep(this,1,length(groups))
}
if(length(pcex)!=length(groups))
{
this<-pcex[1]
pcex<-rep(this,1,length(groups))
}
if(length(pcol)!=length(groups))
{
this<-pcol[1]
pcol<-rep(this,1,length(groups))
}
if(length(lty)!=length(groups))
{
this<-lty[1]
lty<-rep(this,1,length(groups))
}
if(length(lcol)!=length(groups))
{
this<-lcol[1]
lcol<-rep(this,1,length(groups))
}
if(length(tfont)!=length(groups))
{
this<-tfont[1]
tfont<-rep(this,1,length(groups))
}
if(length(tcol)!=length(groups))
{
this<-tcol[1]
tcol<-rep(this,1,length(groups))
}
#
for(j in 1:length(groups))
{
sp<-groups[[j]]
edges<-vector()
dat<-coords[sp,]
names<-data$names[sp]
for(i in 1:dim(alledges)[1])
{
this<-alledges[i,1:2]
if(length(this[this %in% sp])==length(this))
{
edges<-rbind(edges,alledges[i,])
}
}#check out all needed links
for(i in 1:dim(edges)[1])
{
x1<-dat[which(names==edges[i,1]),1]
y1<-dat[which(names==edges[i,1]),2]
x2<-dat[which(names==edges[i,2]),1]
y2<-dat[which(names==edges[i,2]),2]
wei<-edges[i,3]/maxi*5
if(scaled==TRUE)
{
segments(x1,y1,x2,y2,col=lcol[j],lty=lty[j],lwd=wei)
}else
{
segments(x1,y1,x2,y2,col=lcol[j],lty=lty[j])
}
}#i
#add points and texts
points(dat,pch=pch[j],cex=pcex[j],col=pcol[j],bg=bg[j])
if(addlabels==TRUE)
{
text(dat,labels=names,font=tfont[j],col=tcol[j])
}
}#j
}


