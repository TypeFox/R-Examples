VCV.array <- function(phy, dim=2, compact=TRUE){
    
    ## turns a phylogeny into a 3d array similar to a VCV matrix
    ## but keeping each beanch length separate. This is useful for 
    ## handling branch length transformations in functions where
    ## VCVs are used to handle phylogenetic structure
    
    ## rewritten to use new ape, via the clade matrix structure
    
    # an extended version that could replace vcv.phylo.array and vcv.phylo (~ 3x faster than it for one thing)
    
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    
    if(is.null(phy$edge.length))
    	stop("Object \"phy\" is missing edge lengths. Cannot provide VCV matrix.")

    
    if(! dim %in% 2:3) stop("dim must be 2 or 3, for a VCV matrix or array respectively. ")
    
    cm <- clade.matrix(phy)
    cmM <- cm$clade.matrix
    cmE <- cm$edge.length
    cmEM <-  cmM*cmE
    
    if(dim == 2){
        V <- crossprod(cmM, cmEM)
        dimnames(V) <- list(phy$tip.label, phy$tip.label)
    } else {
        if(compact){
            
            nTip <- dim(cmM)[2]
            max.node.depth <- max(colSums(cmM))
    
           V <- array(0, dim=c(nTip, nTip, max.node.depth), dimnames=list(phy$tip.label, phy$tip.label, NULL))
    
            ##  ## must be a way of 'applying' or 'outering' this next bit off the clade matrix
            ##  for(i in 1:nTip){
            ##      for(j in 1:nTip){
            ##      	## get the shared edge lengths and insert into the array
            ##      	shared.edge.lengths <- cmE[as.logical(cmM[,i] * cmM[,j])]
            ##      	V[i,j,seq(along=shared.edge.lengths)] <- shared.edge.lengths
            ##      }
            ##  }
            
            ## trialing a method with only one loop that deals with a matrix slice at a time...
            ## does seem to be ~ 2.3 times faster than the above
            
            for(i in 1:nTip){
                Vslice <- cmEM * cmM[,i]
                Vind <- which(Vslice > 0, arr.ind=TRUE)
                Vval <- Vslice[Vind]
                Vrle <- rle(Vind[,2]) 
                Vind[,1]  <-  unlist(mapply(seq, from=1, Vrle$lengths))
                V[i,,][Vind[,c(2,1)]] <- Vval
            }
            
        } else {
            ## returns a big 3d array showing, for each pair of tips, either 0 (not shared) 
            ## or the appropriate edge length if the node is shared - not good on big trees!
            ## but it is faster than the previous version
            
            ##  multiply each column of the edge length matrix by the clade matrix
            V <- apply(cmEM, 2, function(X) X * cmM)
            dims <- dim(cmM)
        
            ##  gives a (Nnodes * Ntips) by Ntips matrix, which needs reshaping into
            ##  an array of Nnodes by Ntips by Ntips and then rotating to Ntips * Ntips * Nnodes
            V <- array(V, rep(dims, c(1,2)))
            
        V <- aperm(V, c(2,3,1))
        }
    }
    
    class(V) <- "VCV.array"
    return(V)
    
}

## ## time trialing v ape
## sz <- ceiling(2^seq(3, 10, by=0.5))
## tm <- matrix(NA, ncol=3, nrow=length(sz))
## 
## for(t in seq(along=sz)){
##     tree <- rcoal(sz[t])
##     tm[t,1] <- system.time(x1 <- vcv.phylo(tree))[3]
##     tm[t,2] <- system.time(x2 <- VCV.array(tree))[3]
##     tm[t,2] <- system.time(x2 <- VCV.array(tree))[3]
##     tm[t,3] <- system.time(x3 <- VCV.array(tree, dim=3))[3]
##     if(sum(x1 - x2) > 1e-10) stop("disagreeement!")
##     if(sum(x1 - apply(x3, c(1,2), sum)) > 1e-10) stop("disagreeement2!")
## }
## 
## plot(tm[,3] ~ sz, typ="l", log='xy', ylim=range(tm))
## lines(tm[,2] ~ sz, col="red")
## lines(tm[,1] ~ sz, col="blue")

## tm <- as.data.frame(tm)
## names(tm) <- c('ape','caper2d','caper3d')
## tm$size <- sz
## tm <- structure(list(ape = c(0.00100000004749745, 0.00299999990966171, 
## 0.00400000007357448, 0.00800000003073364, 0.0140000000828877, 
## 0.0270000000018626, 0.0489999999990687, 0.0989999999292195, 0.207999999984168, 
## 0.379000000073574, 0.743999999947846, 1.51699999999255, 2.99700000009034, 
## 5.76899999997113, 11.6360000000568), caper2d = c(0.00199999997857958, 
## 0.00300000002607703, 0.00400000007357448, 0.00799999991431832, 
## 0.0109999999403954, 0.0159999999450520, 0.0279999999329448, 0.0470000000204891, 
## 0.0700000000651926, 0.116000000038184, 0.218999999924563, 0.373000000021420, 
## 0.739000000059605, 1.4670000000624, 2.65700000000652), caper3d = c(0.00599999993573874, 
## 0.00899999996181577, 0.0130000000353903, 0.0260000000707805, 
## 0.0419999998994172, 0.0670000000391155, 0.134999999892898, 0.291999999899417, 
## 0.576999999932013, 1.05599999998230, 2.58199999993667, 5.92799999995623, 
## 11.826999999932, 34.3509999999078, 73.8229999999749), size = c(8, 
## 12, 16, 23, 32, 46, 64, 91, 128, 182, 256, 363, 512, 725, 1024
## )), .Names = c("ape", "caper2d", "caper3d", "size"), row.names = c(NA, 
## -15L), class = "data.frame")

