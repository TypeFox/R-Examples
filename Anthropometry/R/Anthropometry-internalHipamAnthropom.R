#All these functions are common to the HIPAM_{MO} and HIPAM_{IMO} algorithms.
initialize.tree <- function(x, maxsplit, orness, type, ah=c(23,28,20,25,25), verbose, ...){
 if(type == "MO"){
  x.ps <- getBestPamsamMO(x, maxsplit, orness, type, ah, verbose, ...)
 }else{
   x.ps <- getBestPamsamIMO(x, maxsplit, orness, type, ah, verbose, ...)
  }
  n.clust <- length(unique(x.ps$clustering))
  list(clustering = x.ps$clustering, asw = x.ps$asw, n.levels = 1, medoids = x.ps$medoids, active = 
         rep(TRUE,n.clust), development = matrix(1:n.clust), n.clust = n.clust, metric = "McCulloch")
}


ext.dist <- function(x, maxsplit, orness, ah, verbose){

 num.variables <- dim(x)[2]
 w <- weightsMixtureUB(orness, num.variables)

 num.personas <- dim(x)[1]
 datosm <- as.matrix(x)
 datost <- aperm(datosm, c(2,1))                     
 dim(datost) <- c(1, num.personas * num.variables)   
 rm(datosm)

 bh <- (apply(as.matrix(log(x)), 2, range)[2,] - apply(as.matrix(log(x)), 2, range)[1,]) / ((maxsplit - 1) * 8)
 if(any(bh==0)){ 
   bh[which(bh == 0)]=0.0000000001
 }
 bl <- -3 * bh
 ah <- ah
 al <- 3 * ah

 DIST <- getDistMatrix(datost, num.personas, num.variables, w, bl, bh, al, ah, verbose) 

  if( is.null(dimnames(x)[1]) ) dimnames(DIST) <- NULL
 
 return(DIST)
}


pamsam <- function(x, k = k, type = type, DIST, maxclust = 50, further.maxclust = 20, exhaustive = FALSE, 
                   max.check = 20, total.check = 100, maxsplit = maxsplit, orness = orness, 
                   ah, verbose){
 
 if(type == "MO"){                 
  #Calculate DIST:
  DIST <- ext.dist(x, maxsplit = k, orness = orness, ah, verbose)
 }else{
   DIST <- DIST
  }

 #Initial clustering:
 init.clust <- initial.clustering(x, k, DIST = DIST, maxclust, exhaustive)
 if(!is.na(k)){
  clustering <- init.clust$clustering
  medoids <- init.clust$medoids
  asw.profile <- init.clust$asw.profile
  sesw.profile <- init.clust$sesw.profile
  asw <- asw.profile
  sesw <- sesw.profile
  furth.info <- NULL
 }else{
   k <- init.clust$num.of.clusters
   asw.profile <- init.clust$asw.profile
   sesw.profile <- init.clust$sesw.profile
   asw <- asw.profile[k]
   sesw <- sesw.profile[k]
   medoids <- init.clust$medoids
   clustering <- init.clust$clustering


   #Further partitioning, ordering and collapsing:
   furth.clust <- further.clustering(x, k, DIST = DIST, medoids, clustering, asw, sesw, max.check, total.check, 
                                     maxclust = further.maxclust, exhaustive, maxsplit = maxsplit, 
                                     orness = orness, ah = ah, verbose)
   medoids <- furth.clust$medoids
   clustering <- furth.clust$clustering
   k <- furth.clust$num.of.clusters
   asw <- furth.clust$asw
   furth.info <- list(k.development = furth.clust$k.info, asw.development = furth.clust$asw.info, 
                      sesw.development = furth.clust$sesw.info, total.checks = furth.clust$total.checks)
  }
  if( !is.null(dimnames(x)[[1]]) ) names(clustering) <- dimnames(x)[[1]]


  #OUTPUT:
  list(medoids = medoids, clustering = clustering, asw = asw, num.of.clusters = k, info = furth.info, 
       profiles = list(asw.profile = asw.profile, sesw.profile = sesw.profile), metric = "McCulloch" )
}




#INITIAL CLUSTERING#
#Search function for finding the best clustering of x using PAM or CLARA.  Uses a grid or exhaustive search up 
#to maxclust unless k is given.
initial.clustering <- function(x, k, DIST, maxclust, exhaustive){

 n.obj <- nrow(x)
 
 DIST.init <- full2dist(DIST)

 x.clust <- cluster.search(x, k, DIST = DIST.init, maxclust, exhaustive)

 list(medoids = x.clust$medoids, clustering = x.clust$clustering,asw.profile = x.clust$asw.profile, 
      sesw.profile = x.clust$sesw.profile, num.of.clusters = x.clust$num.of.clusters)
}

cluster.search <- function(x, k, DIST = NULL, maxclust = 50, exhaustive){
   
 n.obj <- nrow(x)
 #Start:
 if(n.obj <= 2 | !is.na(k)) {  #When k is given or < 3 objects.
  if(n.obj <= 2) {
   x.medoids <- x
   k <- n.obj
   x.clustering <- 1:k
   asw.profile <- NA
   sesw.profile <- NA
   names(asw.profile) <- list( k )
   names(sesw.profile) <- list( k )
  }else{
    x.clust <- pam(DIST, k, diss=TRUE)
    asw.profile <- x.clust$silinfo$avg.width
    sesw.profile <- NA
          
    names(asw.profile) <- list( k )
    names(sesw.profile) <- list( k )
    maxclust <- k
   }
 }else{  #Estimate k using average silhouette width:
   maxclust <- min( maxclust, n.obj - 1 )

   if(maxclust > 12 & exhaustive == FALSE){
    gp <- sqrt(maxclust)
    gpr <- round(gp)
    gl <- floor(gp)
    gst <- floor( ((gpr - 2) + (maxclust - gpr * gl)) / 2 + 2 )
    grid <- seq(gst, by = gpr, length = gl)
    gridout <- c(1, grid, maxclust + 1)
   }else {
     grid <- 2:maxclust
     gpr <- 1
    }
        
   aswmax <- -1
   k <- 1
   asw.profile <- rep(-1,maxclust)
      
   sesw.profile <- rep(NA,maxclust) 
        
    for(i in grid){
     clusti <- pam(DIST, i, diss=TRUE)
     asw.profile[i]  <- clusti$silinfo$avg.width
        
      if(asw.profile[i] > aswmax) {
       aswmax <- asw.profile[i]
       k <- i
       x.clust <- clusti
      }
    }
     
   if(gpr > 1){
    kpos <- match(k,gridout)
    grid2 <- c((gridout[kpos - 1] + 1) : (k - 1), (k + 1) : (gridout[kpos + 1] - 1))

     for(i in grid2) {
      clusti <- pam(DIST, i, diss=TRUE)
      asw.profile[i] <- clusti$silinfo$avg.width
               
       if(asw.profile[i] > aswmax){
        aswmax <- asw.profile[i]
        k <- i
        x.clust <- clusti
       }
      }
    }  
  
  }
 #Check that k is a possible value:
    if(n.obj > 2) {
     if( !is.element(k,2:maxclust) ) {
      stop("Number of clusters found is wrong")
     }
      x.medoids <- x[x.clust$medoids,]
      x.clustering <- x.clust$clustering
    
    }
    list(medoids = x.medoids, clustering = x.clustering,asw.profile = asw.profile, sesw.profile = NA, 
         num.of.clusters = k)
}


##########
#ASW CALC#
#Calculates the average silhouette width of the data given a specific clustering. Note that asw.calc needs 
#to be fed the appropriate medoids.
asw.calc <- function(x, k, clustering, DIST = NULL){ 
  n.obj <- nrow(x)
   if(is.matrix(DIST) && (ncol(DIST)==n.obj & nrow(DIST)==n.obj)){
    col.size <- n.obj
    sil.mem <- clustering
    ind <- rep( clustering, col.size ) + rep( k * (0:(col.size-1)), each=n.obj )
    sumdists <- matrix(tapply(DIST,ind,sum),ncol=col.size)
   }
    
    identif <- matrix(0,k,col.size)
    identif[cbind(sil.mem,1:col.size)] <- 1
    clus.sizes <- tabulate(clustering)
    denominators <- matrix(rep(clus.sizes, col.size), ncol=col.size) - identif
    sil.widths <- sumdists / denominators
    ai <- sil.widths[identif == 1]
    bi <-  apply( matrix(sil.widths[identif == 0] , nrow = k - 1), 2, min)
    si <- (bi - ai) / apply(cbind(bi,ai),1,max)
    si[!is.finite(si)] <- 0
    #[gives 0 as value of a(i) when only 1 in cluster].
    asw <- mean(si)
    sesw <- sqrt( var(si) / col.size )
    list(asw = asw, sesw = sesw)
}


#####################
#FURTHER CLUSTERING#
#Further partitioning and collapsing; uses initial.clustering.
further.clustering <- function(x, k, DIST, medoids, clustering, asw, sesw, max.check, total.check, maxclust, 
                               exhaustive, maxsplit, orness, ah, verbose){
  
 n.obj <- nrow(x)
 itc <- 0
 isc <- 0
 clust.perm <- sample(k)
 k.info <- k
 asw.info <- asw
 sesw.info <- sesw
 prev.groups <- NULL
 prev.groups.index <- NULL
  
  while(isc < k & itc <= total.check){
   itc <- itc + 1
   isc <- isc + 1
   clust.i <- clust.perm[isc]
    if(runif(1) < 0.5) {
     poss.clust <- partition.cluster(x, k, DIST = DIST, clust.i, medoids, clustering, maxclust, exhaustive)
      if(poss.clust$asw <= asw) {
       poss.clust <- collapse.cluster(x, k, DIST = DIST, clust.i, medoids, clustering, asw, sesw, maxclust, 
                                      exhaustive, maxsplit, orness, ah, verbose) 
      }
    }else{
      poss.clust <- collapse.cluster(x, k, DIST = DIST, clust.i, medoids, clustering, asw, sesw, maxclust, 
                                     exhaustive, maxsplit, orness, ah, verbose)
       if(poss.clust$asw <= asw) {
        poss.clust <- partition.cluster(x, k, DIST = DIST, clust.i, medoids, clustering, maxclust, exhaustive)
       }
     }

   if(poss.clust$asw > asw) {
    k <- poss.clust$num.of.clusters
    asw <- poss.clust$asw
    sesw <- poss.clust$sesw
    medoids <- poss.clust$medoids
    clustering <- poss.clust$clustering
    k.info <- c(k.info, k)
    asw.info <- c(asw.info, asw)
    sesw.info <- c(sesw.info, sesw)
    isc <- 0 
    clust.perm <- sample(k)
   }
  }

    list(medoids = medoids, clustering = clustering, num.of.clusters = k, asw = asw, sesw = sesw, k.info = k.info, 
         asw.info = asw.info, sesw.info = sesw.info, total.checks = itc)
}


#PARTITION CLUSTER#
partition.cluster <- function(x, k, DIST = NULL, j, medoids, clustering, maxclust, exhaustive){

 if(length(clustering[clustering==j]) >= 1){
 #Partitioning:
  if(is.null(DIST)) {
   DIST.part <- NULL
  }else{
    DIST.part <- full2dist( DIST[clustering==j,clustering==j] )
   }
    
     clust.obj <- cluster.search(x[clustering==j,], k = NA, DIST = DIST.part, maxclust = maxclust, exhaustive)
     kj <- clust.obj$num.of.clusters
     medoids.j <- clust.obj$medoids
    
  #Ordering:
  clustering[clustering > j] <- clustering[clustering>j] + kj - 1
  clustering[clustering == j] <- j + clust.obj$clustering - 1
  if(j == 1)        medoids <- rbind(medoids.j, medoids[(j+1):k,])
  if(j > 1 & j < k) medoids <- rbind(medoids[1:(j-1),], medoids.j, medoids[(j+1):k,])
  if(j == k)        medoids <- rbind(medoids[1:(j-1),], medoids.j)

  k <- k + kj - 1 
  asw.obj <- asw.calc(x, k, clustering, DIST)
  asw <- asw.obj$asw
  sesw <- asw.obj$sesw
 }
 
 list(medoids = medoids, clustering = clustering,asw = asw, sesw = sesw, num.of.clusters = k)
}


##################
#COLLAPSE CLUSTER#
collapse.cluster <- function(x, k, DIST, j, medoids, clustering, asw, sesw, maxclust, exhaustive, maxsplit, 
                             orness, ah,verbose){

 if(k > 2){ 
  DIST.medoids <- ext.dist(medoids, maxsplit, orness = orness, ah, verbose) 
  cluster.weights <- 1/(DIST.medoids[j,-j] + .001)
  cluster.weights <- cluster.weights/sum(cluster.weights)
  chosen.order <- sample((1:k)[-j], prob = cluster.weights)
 } 
  n.check <- min(k-1,3)
  ip <- 0
   while(k > 2 & ip < n.check){ 
    ip <- ip + 1
    change1 <- j
    change2 <- chosen.order[ip]
    coll.step <- collapse.step(x, k, DIST, change1, change2, medoids, clustering, asw, sesw)
     if(coll.step$asw > asw) {
      medoids <- coll.step$medoids
      clustering <- coll.step$clustering
      asw <- coll.step$asw
      sesw <- coll.step$sesw
      k <- coll.step$num.of.clusters
      ip <- n.check + 1
     }
   }
  
  list(medoids = medoids, clustering = clustering,num.of.clusters = k, asw = asw, sesw = sesw)
}


#COLLAPSE STEP#
collapse.step <- function(x, k, DIST, change1, change2, medoids, clustering, asw, sesw){
 ch.low <- min(change1,change2)
 ch.high <- max(change1,change2)
  if(ch.high != ch.low){
   poss.clustering <- clustering
   poss.clustering[poss.clustering == ch.high] <- ch.low
   poss.clustering[poss.clustering > ch.high] <- poss.clustering[poss.clustering > ch.high] - 1
   poss.medoids <- medoids
   poss.medoids <- medoids[-ch.high,]
   x.mm <- x[poss.clustering==ch.low,]
      DIST.mm <- DIST[poss.clustering == ch.low, poss.clustering == ch.low]
      pos.mm <- pam(full2dist(DIST.mm), k=1)$medoids
      poss.medoids[ch.low,] <- x.mm[pos.mm,]
       if(!is.null(dimnames(medoids))) dimnames(medoids)[[1]][ch.low] <- pos.mm
       k.poss <- k - 1 
       asw.poss.obj <- asw.calc(x, k=k.poss, poss.clustering, DIST) 
       asw.poss <- asw.poss.obj$asw
       sesw.poss <- asw.poss.obj$sesw
       if(asw.poss >= asw){
        asw <- asw.poss
        sesw <- sesw.poss
        clustering <- poss.clustering
        k <- k.poss
        medoids <- poss.medoids
       }
 }else{
  }
    
  list(medoids = medoids, clustering = clustering,num.of.clusters = k, asw = asw, sesw = sesw)
}

hipam.local <- function(tree, x, asw.tol, maxsplit, local.const, orness, type, ah, verbose,...){
 #Iterative procedure to keep splitting branches as long as there are active branches.
 active <- tree$active
 while(sum(active) > 0){
  choose.among <- (1:length(active))[active]
   if(length(choose.among) == 1){
    whch.brnch <- choose.among
   }else {
     whch.brnch <- sample(choose.among, 1)
    }

   if(type == "MO"){
    proposal <- checkBranchLocalMO(tree, x, whch.brnch, maxsplit, asw.tol, local.const, 
                                   orness, type, ah, verbose, ...)
   }else{
     proposal <- checkBranchLocalIMO(tree, x, whch.brnch, maxsplit, asw.tol, local.const, 
                                     orness, type, ah, verbose, ...)
    }
 
    proposed.tree <- proposal$tree
    if(proposal$reject){
     tree$active[whch.brnch] <- FALSE
    }else{
       tree <- proposed.tree
      }
      active <- tree$active
 }

    tree
}


update.tree.local <- function(object, xi.ps, which.x, i){
 tree <- object
 #Number of clusters before the update:
 maxclust <- max(tree$clustering)
 
 #Update clustering:
 tree$clustering[which.x] <- maxclust + xi.ps$clustering
 
 #Update number of clusters:
 new.clust <- xi.ps$num.of.clusters
 tree$n.clust <- tree$n.clust + new.clust - 1
 clust <- tree$clustering
 
 #Update medoids:
 tree$medoids <- rbind(tree$medoids, xi.ps$medoids)
 
 #Update active status:
 tree$active[i] <- FALSE
 new.activities <- rep(TRUE, new.clust)
 
  for (j in 1:new.clust){
   if (sum(xi.ps$clustering==j) <= 2){
    new.activities[j] <- FALSE
   }
  }
  
  tree$active <- c(tree$active,new.activities)
  
  #Update the development matrix:
  dmat<- tree$development
  nr <- nrow(dmat)
  nc <- ncol(dmat)
  whch.row <- sum((1:nr) * apply((dmat==i), 1, sum, na.rm = TRUE))
  whch.col <- sum((1:nc) * apply((dmat==i), 2, sum, na.rm = TRUE))
  dmat<-dmat[c(1:whch.row,rep(whch.row,new.clust-2),whch.row:nr),]
  
  if (nc==whch.col){
   dmat <- cbind(dmat,NA)
   tree$n.levels <- tree$n.levels + 1
  }
  
  for (j in 1:new.clust){
   dmat[whch.row + j - 1, whch.col + 1] <- maxclust + j
  }

 tree$development <- dmat
 tree
}


full2dist <- function(fullm, method = "unspecified"){
    distm <- fullm[lower.tri(fullm)] #?lower.tri: Returns a matrix of logicals the same size of a given 
    #matrix with entries. 
    #TRUE in the lower or upper triangle. 
    attr(distm , "Size") <- nrow(fullm)
    attr(distm , "Labels") <- dimnames(fullm)[[1]]
    attr(distm , "Diag") <- FALSE
    attr(distm , "UPPER") <- FALSE
    attr(distm, "method") <- method
    attr(distm , "class") <- "dist"
    attr(distm , "call") <- match.call()
    class(distm) <- "dist"
    return(distm)
}
