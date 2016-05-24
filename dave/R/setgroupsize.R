# FUNCTION SETGROUPSIZE (internal use only)
# =========================================
setgroupsize<- function(vec){
# does what table() does, but without sorting the factors
# 29. 5. 2004 - na's now counted as groups
vec<- as.integer(vec)
redgr<- 0
vec[is.na (vec)] <- 0
ngroups<- length(table(vec))
grouplabs<- rep(0,ngroups)
groupcounts<- rep(0,ngroups)
vlen<- length(vec)
grouplabs[1]<- vec[1]
groupcounts[1]<- 1
labind<- 1
for (i in 2:vlen){
if(vec[i] == 0) redgr<- 1
if(vec[i] == grouplabs[labind]) groupcounts[labind]<- groupcounts[labind]+1
if(vec[i] != grouplabs[labind]){
labind<- labind+1 
groupcounts[labind]<- groupcounts[labind]+1
grouplabs[labind]<- vec[i]
}
}
outsetgroupsize<- list(ngroups=ngroups,grouplabs=grouplabs,groupcounts=groupcounts,vlen=vlen,vec=vec)
}
