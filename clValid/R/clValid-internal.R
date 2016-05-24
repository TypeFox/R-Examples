####################################################################
## Internal functions for clValid
####################################################################


####################################################################
## vClusters() function
####################################################################
## Arguments
####################################################################
## mat - data matrix
## clMethod - clustering method
## nClust - minimun number of clusters to evaluate
## nclustMax - maximun number of clusters to evaluate
## ... - arguments to pass to clustering functions
####################################################################
## Value
####################################################################
## List with:
## clusterObj - clustering results
## measures - validation measures
####################################################################


vClusters <- function(mat,clMethod,nClust,nclustMax, validation,
                      Dist, method, metric, annotation, GOcategory,
                      goTermFreq, neighbSize, dropEvidence, verbose, ... ) {
  
  measNames <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                 if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                 if("biological"%in%validation) c("BHI","BSI"))
  measures <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  rownames(measures) <- measNames
  colnames(measures) <- nClust

  switch(clMethod,
         hierarchical = {
           clusterObj <- hclust(Dist,method)
         },
         diana = {
           clusterObj <- diana(Dist, ...)
         },
         kmeans = {
           clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust
           clusterObjInit <- hclust(Dist,method)
         },
         agnes = {
           clusterObj <- agnes(Dist, method=method, ...)
         },
         ## otherwise - sota, fanny, som, model, pam, clara
         { clusterObj <- vector("list",length=length(nClust))
           names(clusterObj) <- nClust })

  ind <- 1
  for (nc in nClust) {
    switch(clMethod,
           kmeans = {
             initial <- tapply(mat, list(rep(cutree(clusterObjInit,nc),ncol(mat)),col(mat)),
                               function(x) mean(x, na.rm=TRUE))
             if(length(dup <- which(duplicated(initial)))>0) {
               for(dupi in dup) 
                 initial[dupi,] <- initial[dupi,] + jitter(initial[dupi,])
             }
             dimnames(initial) <- list(NULL,dimnames(mat)[[2]])
             clusterObj[[ind]] <- kmeans(mat,initial,...)
             cluster <- clusterObj[[ind]]$cluster
           },
           fanny = {
             clusterObj[[ind]] <- fanny(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           model = {
             clusterObj[[ind]] <- Mclust(mat,nc, ...)
             cluster <- clusterObj[[ind]]$classification
           },
           som = {
             clusterObj[[ind]] <- som(mat, grid=somgrid(1,nc), ...)
             cluster <- clusterObj[[ind]]$unit.classif
           },
           pam = {
             clusterObj[[ind]] <- pam(Dist, nc, ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           clara = {
             clusterObj[[ind]] <- clara(mat, nc, metric=ifelse(metric=="correlation","euclidean",metric), ...)
             cluster <- clusterObj[[ind]]$clustering
           },
           sota = {
             clusterObj[[ind]] <- sota(mat,nc-1)
             cluster <- clusterObj[[ind]]$clust

           },
           ## otherwise - hierarchical, diana, agnes
           {cluster <- cutree(clusterObj,nc)})

    if(length(table(cluster))!=nc) {
      warning(paste(clMethod, "unable to find",nc,"clusters, returning NA for these validation measures"))
      measures[,ind] <- NA
      ind <- ind+1
      next()
    }
    
    ## internal validation measures
    if ("internal"%in%validation) {
      measures["Dunn",ind] <- dunn(Dist ,cluster)
      measures["Silhouette",ind] <- mean(silhouette(cluster, dmatrix=as.matrix(Dist))[,3])
      measures["Connectivity",ind] <- connectivity(Dist ,cluster, neighbSize=neighbSize)
      if(verbose) print(paste("Finished internal validation,", clMethod, nc, "clusters"))
    }
    
    if("biological"%in%validation) {
      measures["BHI",ind] <- BHI(cluster,annotation=annotation, names=rownames(mat),
                                 category=GOcategory, dropEvidence=dropEvidence)
      if(verbose & "biological"%in%validation)
        print(paste("Finished BHI,", clMethod, nc, "clusters"))      
    }

    ## stability validation measures
    if ("stability"%in%validation | "biological"%in%validation) {
      co.del <- 0 ## for use in verbose printing of progress
      for (del in 1:ncol(mat)) {
        matDel <- mat[,-del]               ## matDel <- as.matrix(matDel)
        if(metric=="correlation") {
          DistDel <- as.dist(1-cor(t(matDel), use="pairwise.complete.obs"))
        } else {
          DistDel <- dist(matDel,method=metric)
        }
        switch(clMethod,
               hierarchical = clusterObjDel <- hclust(DistDel,method),
               kmeans = clusterObjInitDel <- hclust(DistDel,method),
               diana = clusterObjDel <- diana(DistDel, ...),
               agnes = clusterObjDel <- agnes(DistDel, method=method, ...),
               clara = clusterObjDel <- clara(matDel,nc,metric=ifelse(metric=="correlation","euclidean",metric), ...))


        switch(clMethod,
               kmeans = {
                 initialDel <- tapply(matDel, list(rep(cutree(clusterObjInitDel,nc),
                                                       ncol(matDel)), col(matDel)),
                                      function(x) mean(x, na.rm=TRUE))
                 if(length(dup <- which(duplicated(initialDel)))>0) {
                   for(dupi in dup) 
                     initialDel[dupi,] <- initialDel[dupi,] + jitter(initialDel[dupi,])
                 }
                 dimnames(initialDel) <- list(NULL,dimnames(matDel)[[2]])
                 kmdel <- kmeans(matDel,initialDel, ...)
                 clusterDel <- kmdel$cluster
               },
               fanny = {
                 hfdel <- fanny(DistDel, nc, ...)
                 clusterDel <- hfdel$clustering
               },
               model = {
                 clusterDel <- Mclust(matDel,nc, ...)$classification
               },
               som = {
                 hsdel <- try(som(matDel, grid=somgrid(1,nc), ...))
                 clusterDel <- hsdel$unit.classif
               },
               pam = {
                 clusterDel <- pam(DistDel, nc, cluster.only=TRUE, ...)
               },
               clara = {
                 clusterDel <- clusterObjDel$clustering
               },
               sota = {
                 clusterDel <- sota(matDel,nc-1)$clust
               },
               ## otherwise - hierarchical, diana, agnes
               {clusterDel <- cutree(clusterObjDel,nc)})

        if("stability"%in%validation) {
          stabmeas <- stability(mat, Dist, del, cluster, clusterDel)
          measures["APN",ind] <- measures["APN",ind] + stabmeas["APN"]
          measures["AD",ind]  <- measures["AD",ind]  + stabmeas["AD"]
          measures["ADM",ind] <- measures["ADM",ind] + stabmeas["ADM"]
          measures["FOM",ind] <- measures["FOM",ind] + stabmeas["FOM"]
        }
        if("biological"%in%validation) {
          tmp <- BSI(cluster,clusterDel,annotation=annotation,
                     names=rownames(mat), category=GOcategory, goTermFreq=goTermFreq,
                     dropEvidence=dropEvidence)
          measures["BSI",ind] <- measures["BSI",ind] + tmp
        }
        ## VERBOSE printing
        if (del/ncol(mat) > 0.25 & co.del==0) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 25% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 25% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
        else if (del/ncol(mat) > 0.50 & co.del==1) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 50% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 50% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
        else if (del/ncol(mat) > 0.75 & co.del==2) 
          {
            if(verbose & "stability"%in%validation) 
              print(paste("Stability validation 75% finished,", clMethod, nc, "clusters"))
            if(verbose & "biological"%in%validation) 
              print(paste("BSI 75% finished,", clMethod, nc, "clusters"))            
            co.del <- co.del+1
          }
      } #END OF del LOOP
      if(verbose & "stability"%in%validation)
        print(paste("Finished stability validation,", clMethod, nc, "clusters"))
      if(verbose & "biological"%in%validation)
      print(paste("Finished BSI,", clMethod, nc, "clusters"))
    } #END of STABILITY measures
    ind <- ind+1  #ind tracks number clusters
    ## if(verbose) print(paste("Finished with", nc, "clusters"))
  } #END OF NC LOOP
  
  if ("stability"%in%validation) {
    measures["APN",] <- measures["APN",]/ncol(mat)
    measures["AD",] <-  measures["AD",]/ncol(mat)
    measures["ADM",] <- measures["ADM",]/ncol(mat)
    measures["FOM",] <- measures["FOM",]/ncol(mat)  ## little different from Yeung paper (doesn't do this)
  }
  if ("biological"%in%validation) {
    measures["BSI",] <- measures["BSI",]/ncol(mat)
  }
  
  list(clusterObj=clusterObj, measures=measures)
}

