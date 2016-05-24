
#
#  This code was used to create structures defining the polyeders/
#  triangulations used in 3D visualization
#
phi <- (1+sqrt(5))/2
icosa0 <- list(vertices=matrix(c(0,1,phi,
                             0,-1,phi,
                             0,1,-phi,
                             0,-1,-phi,
                             1,phi,0,
                             1,-phi,0,
                             -1,phi,0,
                             -1,-phi,0,
                             phi,0,1,
                             phi,0,-1,
                             -phi,0,1,
                             -phi,0,-1),3,12)/sqrt((5+sqrt(5))/2),
                    indices=c(1,2,9,
                              1,2,11,
                              1,5,7,
                              1,5,9,
                              1,7,11,
                              2,6,8,
                              2,6,9,
                              2,8,11,
                              3,4,10,
                              3,4,12,
                              3,5,7,
                              3,5,10,
                              3,7,12,
                              4,6,8,
                              4,6,10,
                              4,8,12,
                              5,9,10,
                              6,9,10,
                              7,11,12,
                              8,11,12),
                    edges=matrix(c(1,2,
                                   1,5,
                                   1,7,
                                   1,9,
                                   1,11,
                                   2,6,
                                   2,8,
                                   2,9,
                                   2,11,
                                   3,4,
                                   3,5,
                                   3,7,
                                   3,10,
                                   3,12,
                                   4,6,
                                   4,8,
                                   4,10,
                                   4,12,
                                   5,7,
                                   5,9,
                                   5,10,
                                   6,8,
                                   6,9,
                                   6,10,
                                   7,11,
                                   7,12,
                                   8,11,
                                   8,12,
                                   9,10,
                                   11,12),2,30),
                              nv=12,
                              ni=20,
                              ne=30)
refine.polyeder <- function(polyeder){
norm <- function(x) sqrt(sum(x^2))
ni <- polyeder$ni
nv <- polyeder$nv
ne <- polyeder$ne
vertices <- polyeder$vertices
indices <- matrix(polyeder$indices,3,ni)
edges <- polyeder$edges
nnv <- nv+ne
nni <- 4*ni
nne <- nnv+nni-2# ne+ne+3*ni
nvertices <- matrix(0,3,nnv)
nindices <- array(0,c(3,4*ni))
nedges <- matrix(0,2,nne)
nvertices[,1:nv] <- vertices
vorigin <- matrix(0,2,ne)
for(i in 1:ne){
   vert <- apply(vertices[,edges[,i]],1,mean)
   e1 <- edges[,i]
   v1 <- i+nv
   nedges[,2*i-1] <- c(e1[1],i+nv)
   nedges[,2*i] <- c(e1[2],i+nv)
   vorigin[,i] <- e1
   nvertices[,nv+i] <- vert/norm(vert)
}
ne2 <- 2*ne
ne22 <- ne2-2
ne21 <- ne2-1
for(i in 1:ni){
   i4 <- 4*i
   i3 <- 3*i
   ind0 <- indices[,i]
   e1 <- ind0[1:2]
   e2 <- ind0[c(1,3)]
   e3 <- ind0[2:3]
   v1 <- (1:ne)[apply((vorigin-e1)==0,2,all)]+nv
   v2 <- (1:ne)[apply((vorigin-e2)==0,2,all)]+nv
   v3 <- (1:ne)[apply((vorigin-e3)==0,2,all)]+nv
   vsort <- sort(c(v1,v2,v3))
   v1 <- vsort[1]
   v2 <- vsort[2]
   v3 <- vsort[3]
   nindices[,i4-3] <- c(ind0[1],v1,v2)
   nindices[,i4-2] <- c(ind0[2],v1,v3)
   nindices[,i4-1] <- c(ind0[3],v2,v3)
   nindices[,i4] <- c(v1,v2,v3)
   nedges[,ne22+i3] <- c(v1,v2)
   nedges[,ne21+i3] <- c(v1,v3)
   nedges[,ne2+i3] <- c(v2,v3)
}
list(vertices=nvertices,indices=as.vector(nindices),edges=nedges,nv=nnv,ni=nni,ne=nne)
}  
icosa1 <- refine.polyeder(icosa0)
icosa2 <- refine.polyeder(icosa1)
icosa3 <- refine.polyeder(icosa2)
icosa4 <- refine.polyeder(icosa3)
icosa5 <- refine.polyeder(icosa4)
save(icosa0,icosa1,icosa2,icosa3,icosa4,icosa5,file="polyeders.rda",compress = TRUE)
#
#  the file "polyeders.rda" needs to be copied to dti/data/
#