
# Conditionally Independent Dyadic Models for social and other networks
# CID Team, CMU Statistics
# ACT edited this Oct 13 2013

#library(Rcpp); library(mvtnorm); library(msm); sourceCpp("cid.cpp")

#class members: initialize, reinitialize,  value, value.ext,  generate, pieces, log.likelihood, draw,  gibbs.full, gibbs.value
#potential members: random.start

#redo the C++ function listing here, just in case. Build only.
#library(Rcpp); compileAttributes ("../../CIDnetworks")

rdirichlet.one <- function(one) {   #blocks-by-nodes
  r1 <- rgamma(length(one), one)
  r1/sum(r1)
}

rdirichlet.block <- function(node.block) {   #blocks-by-nodes
  r1 <- matrix(rgamma(length(c(node.block)), node.block), nrow=nrow(node.block))
  r1 <- t(r1)*apply(r1,2,max)
  t(r1/rowSums(r1))
}

#Log of the inverse gamma density.
my.dinvgamma <- function (x, a, b) -a*log(b) - lgamma(a) - (a+1)*log(x) - b/x



non.diag <- function(nn) as.logical(1-diag(nn))

l.diag <- function(nn) {
  out <- NULL
  for (kk in 1:(nn-1)) out <- c(out,(kk-1)*nn+(kk+1):nn)
  return(out)
}

u.diag <- function(nn) {
  out <- NULL
  for (kk in 1:(nn-1)) out <- c(out,seq(kk+kk*nn,nn*nn,by=nn))
  return(out)
}

#make ordinal data from basic
ordinal.maker <- function (vec, cuts=quantile(vec, c(0.25, 0.5, 0.75))) {
  #vec=0:20; cuts=cutoffs
  apply(1*outer (vec, cuts, ">="), 1, sum)
}


my.pmvnorm <- function (lower, #n-by-2
                        upper, #n-by-2
                        meanval, #n-by-2
                        sigma, #1
                        rho,
                        log=TRUE) { #1

  ## lower=matrix(c(-2:2/2, -2:2/2), ncol=2); upper=1+lower; meanval=0; sigma=1; rho=0.2
  ## lower=matrix(rep(-Inf, 10), ncol=2); upper=matrix(seq(-5,5,length=10), ncol=2); meanval=0; sigma=1; rho=0.2
  ## lower=matrix(seq(-5,5,length=10), ncol=2); upper=matrix(rep(Inf,10), ncol=2); meanval=0; sigma=1; rho=0.2
  
  new.lower <- (lower - meanval)/sigma
  new.upper <- (upper - meanval)/sigma

  new.lower[new.lower>50] <- 50
  new.lower[new.lower< -50] <- -50
  
  new.upper[new.upper>50] <- 50
  new.upper[new.upper< -50] <- -50
  
  xx <- c(new.upper[,1], new.lower[,1], new.upper[,1], new.lower[,1])
  yy <- c(new.upper[,2], new.upper[,2], new.lower[,2], new.lower[,2])

  r1 <- matrix(pbivnorm(xx,yy,rho), ncol=4)  #library(pbivnorm)
  if (any(is.na(r1))) warning ("Missing values returned from bivariate normal CDF.")
  r1[,2:3] <- -r1[,2:3]
  
  output <- apply(r1, 1, sum)
  if (log) return(log(output)) else return(output)
}



#Link to C++ function.

make.arc.list <- function (nn=10) makeArcList(nn)
make.edge.list <- function (nn=10) makeEdgeList(nn)
#includes self-edges.
make.edge.list.selfies <- function (nn=10) makeEdgeListSelfies(nn)

#Link to C++ function.
edge.list.distance <- function(latent.space.pos, edge.list) {
  out <- eldc(latent.space.pos, edge.list)
  return(out)
}
cosine.closeness <- function(latent.space.pos, edge.list) {
  out <- cosineClosenessC(latent.space.pos, edge.list)
  return(out)
}

#row assignments within edge.list. Speeds up later computation on node level.
row.list.maker <- function (edge.list, nodes=max(c(edge.list)))
  lapply(1:nodes, function(rr) sort(c(which (edge.list[,1] == rr),
                                      which (edge.list[,2] == rr))))

#Remove first "burnin" columns from Gibbs objects.
remove.burnin <- function(gibbs.out, burnin=100) lapply (gibbs.out, function(gg) rbind(rbind(gg)[,-(1:burnin)]))


#X'X for sparse ones matrices: coincident rows.
coincident.rows <- function (rows.list) {
  out <- sapply (1:length(rows.list), function(rr) 
                 sapply (1:length(rows.list), function(ss) 
                         length(intersect(rows.list[[rr]], rows.list[[ss]]))
                         )
                 )
  return(out)
}

xtx <- function (edge.list, nodes) coincidence (edge.list, as.integer(nodes))
xtyc <- function (edge.list, outcome, nodes) xty (edge.list, outcome, as.integer(nodes))

unittest.xtx <- function () {
  nodes=10; edge.list <- make.edge.list(nodes)
  c1 <- xtx(edge.list, nodes)
  c2 <- xtyc(edge.list, rnorm(dim(edge.list)[1]), nodes)

}

# For post-processing outputs, matrices may be preferable to lists.
list.output.to.matrices.simple <- function (gibbs.out) {
  elements <- lapply (1:length(gibbs.out[[1]]), function(el)
                      sapply(gibbs.out, function(it) it[[el]]))
  names(elements) <- names(gibbs.out[[1]])
  elements
}


#Now works for nested lists, like we would find in the CID class.
list.output.to.matrices <- function (gibbs.out) {
  #gibbs.out = test.lsm.g.gibbs
  elements <- lapply (1:length(gibbs.out[[1]]), function(el)
                      if (!grepl("out",class(gibbs.out[[1]][[el]]))) {
                        sapply(gibbs.out, function(it) it[[el]])
                      } else {
                        subval1 <- lapply(1:length(gibbs.out[[1]][[el]]), function (subel) sapply(gibbs.out, function(it) it[[el]][[subel]]))
                        names(subval1) <- paste(class(gibbs.out[[1]][[el]]),names(gibbs.out[[1]][[el]]),sep=".")
                        subval1
                      })
  
  names(elements) <- names(gibbs.out[[1]])
  count <- 0
  elements2 <- list()
  name.vector <- c()
  for (ee in 1:length(gibbs.out[[1]])) 
    if (class(elements[[ee]]) != "list") {count <- count+1; elements2[[count]] <- elements[[ee]]; name.vector[count] <- names(elements)[ee]} else for (kk in 1: length(elements[[ee]])) {count <- count+1; elements2[[count]] <- elements[[ee]][[kk]]; name.vector[count] <- names(elements[[ee]])[kk]}
  names(elements2) <- name.vector
  elements2
}


#postprocess.latent.positions is now more general than one particular family.

#Rotate LSM positions so that nth point lies in the 1:nth dimensional subspace.
postprocess.latent.positions <- function (latent.space.coords, latent.space.target,recenter=TRUE) {
  #latent.space.coords=latent.space.pos

  #recenter data.
  if (recenter){
    latent.space.coords <- t(t(latent.space.coords)-apply(latent.space.coords,2,mean))
    latent.space.target <- t(t(latent.space.target)-apply(latent.space.target,2,mean))
	}
  
  projection = t(latent.space.target)%*% latent.space.coords;
  ssZ = svd(projection)
  transformation = ssZ$v%*%t(ssZ$u)
  latent.space.coords = latent.space.coords%*%transformation
  
  
  #rotate the nth position to the nth axis in plane for each.

#   rot.mat <- function (lsp.row, pinned.d=1) {
#     #lsp.row=latent.space.pos[1,]; pinned.d=1
#     rot.out <- diag(1, length(lsp.row))
#     for (dd in (length(lsp.row):(pinned.d+1))) {
#       d1 <- dd-1; d2 <- dd
#       rot <- diag(1, length(lsp.row))
#       ss <- sqrt(sum(lsp.row[c(d1,d2)]^2))
#       rot[d1,d1] <- rot[d2,d2] <- lsp.row[d1]/ss
#       rot[d1,d2] <- -lsp.row[d2]/ss
#       rot[d2,d1] <- lsp.row[d2]/ss
#       lsp.row <- lsp.row%*%rot
#       rot.out <- rot.out%*%rot
#     }
#     return(rot.out)
#   }
# 
#   if (dim(latent.space.coords)[2] > 1) {
#     for (kk in 1:(dim(latent.space.coords)[2]-1)) {   #rotate into primary (hyper)plane.
#       latent.space.coords <-
#         latent.space.coords%*%rot.mat(latent.space.coords[kk,], kk)
#     }
#     for (kk in 2:dim(latent.space.coords)[2]) {
#       latent.space.coords[,kk] <-
#         latent.space.coords[,kk]*(1*(latent.space.coords[kk,kk] >= 0) - 1*(latent.space.coords[kk,kk] < 0))
#     }
#   }
  
  return(latent.space.coords)
}

