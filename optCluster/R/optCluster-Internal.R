########## optCluter Internal Function ##############

## mbCountVal
## Determine validation scores for count based clustering algorithms
mbCountVal <- function(obj, Normalizer, uniqueTreatment, validation, measNames, nClust, metric, neighbSize, clVerbose, 
 					mbCountMethods, iter.max, TMP, annotation, GOcategory, goTermFreq, dropEvidence){                
		if(!requireNamespace("MBCluster.Seq")) {
      		stop("package 'MBCluster.Seq' required for clustering count data")
    	}

		if(metric=="correlation"){
    		Dist <- as.dist(1-cor(t(obj), use="pairwise.complete.obs"))
    		}  else {
  			Dist <- dist(obj,method=metric)}
  		
  		if("biological" %in% validation){	
  			if(is.character(annotation) & length(grep(".db", annotation))==0) {
    			annotation <- paste(annotation, ".db", sep="")
  			}
  			## Convert annotation list to annotation table
  			if (is.list(annotation)) {
    			if(is.null(rownames(obj))) {
      				stop("rownames of data must be present to specify biological annotation from file")
    			}
    			annotation <- annotationListToMatrix(annotation, genenames=rownames(obj))
  			}

  			if (is.null(annotation)) {
    			stop("annotation must be specified in order to use biological validation")
  			}
  			if (is.character(annotation)) {
    			if(!requireNamespace("Biobase") | !requireNamespace("GO.db") | !requireNamespace("annotate")) {
      				stop("packages 'Biobase', 'GO.db', and 'annotate' required for 2nd type of biological validation \n
					these can be downloaded from Bioconductor (www.bioconductor.org)")
    			}
  			}
  			if(is.null(rownames(obj))) {
    			stop("rownames of data must be present to use biological validation")
  			}			
		}
		## Compute Validation Measures for RNA-Seq clustering
		mbCountData <- MBCluster.Seq::RNASeq.Data(obj, Normalizer = Normalizer, Treatment = uniqueTreatment, GeneID <- rownames(obj))
  		allMeasures <- array(dim=c(length(measNames),length(nClust),length(mbCountMethods)))
  		dimnames(allMeasures) <- list(measNames,nClust,mbCountMethods)
  		colnames(allMeasures) <- nClust

		allClusterObj <- vector("list",length=length(mbCountMethods))
		names(allClusterObj) <- mbCountMethods 
        
        for(i in 1:length(mbCountMethods)){
        	
        	switch(mbCountMethods[i],
        		em.nbinom = {
        			method = "EM"
        			countModel = "nbinom"
        		},
        		da.nbinom = {
        			method = "DA"
        			countModel = "nbinom"
        		},
        		sa.nbinom = {
        			method = "SA"
        			countModel = "nbinom"
        		},
        		em.poisson = {
        			method = "EM"
        			countModel = "poisson"
        		},
        		da.poisson = {
        			method = "DA"
        			countModel = "poisson"
        		},
        		sa.poisson = {
        			method = "SA"
        			countModel = "poisson"
        		}
        		)
        	
        	MBmeasures <- matrix(0,nrow=length(measNames),ncol=length(nClust))
  			rownames(MBmeasures) <- measNames
  			colnames(MBmeasures) <- nClust

			MBclusterObj <- vector("list",length=length(nClust))
			names(MBclusterObj) <- nClust         	           
			ind <- 1
  			for (nc in nClust) {
				centers <- MBCluster.Seq::KmeansPlus.RNASeq(mbCountData, nK = nc)$centers
				noPrint <- capture.output(rnaClust <- MBCluster.Seq::Cluster.RNASeq(mbCountData, model = countModel, centers = centers, method = method,
				iter.max = iter.max, TMP = TMP))
				MBclusterObj[[ind]] <- rnaClust
        		MBcluster <- MBclusterObj[[ind]]$cluster        
        		## Avoid errors in rank aggregation
    			repCount <- 0
    			while(length(table(MBcluster))!=nc) {
					centers <- MBCluster.Seq::KmeansPlus.RNASeq(mbCountData, nK = nc)$centers
					noPrint <- capture.output(rnaClust <- MBCluster.Seq::Cluster.RNASeq(mbCountData, model = countModel, centers = centers, method = method,
						iter.max = iter.max, TMP = TMP))
					MBclusterObj[[ind]] <- rnaClust
        			MBcluster <- MBclusterObj[[ind]]$cluster
					repCount <- repCount + 1
					if(repCount == 5){					
						stop(mbCountMethods[i]," unable to find ",nc," clusters, rank aggregation cannot be performed")
					}
      			}
      			names(MBcluster) <- rownames(obj)      			

    			## internal validation measures
    			if ("internal"%in%validation) {
      				MBmeasures["Dunn",ind] <- dunn(Dist ,MBcluster)
      				MBmeasures["Silhouette",ind] <- mean(silhouette(MBcluster, dmatrix=as.matrix(Dist))[,3])
      				MBmeasures["Connectivity",ind] <- connectivity(Dist ,MBcluster, neighbSize=neighbSize)
      
      			if(clVerbose) print(paste("Finished internal validation", mbCountMethods[i], nc, "clusters"))
    			}
    
    			if("biological"%in%validation) {
      			MBmeasures["BHI",ind] <- BHI(MBcluster,annotation=annotation, names=rownames(obj),
                                 		category=GOcategory, dropEvidence=dropEvidence)
      
        		if(clVerbose & "biological"%in%validation)
        		print(paste("Finished BHI", mbCountMethods[i], nc, "clusters"))      
    			}

    			## stability validation measures
    			if ("stability"%in%validation | "biological"%in%validation) {
      				co.del <- 0 ## for use in verbose printing of progress
      				for (del in 1:ncol(obj)) {
        				objDel <- obj[,-del]
        				if(metric=="correlation") {
          					DistDel <- as.dist(1-cor(t(objDel), use="pairwise.complete.obs"))
        				} else {
          					DistDel <- dist(objDel,method=metric)
        				}
						## identify remaining columns
						uniqueTreatmentDel <- uniqueTreatment[-del]
						if(is.vector(Normalizer)){
							normalizerDel = Normalizer[-del]
						} else if(is.matrix(Normalizer)){
							normalizerDel = Normalizer[,-del]
						} else {
							normalizerDel = NULL
						}
					
						DataDel <- MBCluster.Seq::RNASeq.Data(objDel, Normalizer = normalizerDel, Treatment = uniqueTreatmentDel, GeneID <- rownames(objDel))		
						centersDel <- MBCluster.Seq::KmeansPlus.RNASeq(DataDel, nK = nc)$centers
						noPrintDel <- capture.output(clusterDel <- MBCluster.Seq::Cluster.RNASeq(DataDel, model = countModel, centers = centersDel, 
						method = method, iter.max = iter.max, TMP = TMP)$cluster)
						## Avoid errors in BSI validation
						repCountDel <- 0
    					while(length(table(clusterDel))!=nc) {
							centersDel <- MBCluster.Seq::KmeansPlus.RNASeq(DataDel, nK = nc)$centers
							noPrintDel <- capture.output(clusterDel <- MBCluster.Seq::Cluster.RNASeq(DataDel, model = countModel, centers = centersDel, 
							method = method, iter.max = iter.max, TMP = TMP)$cluster)
							repCountDel <- repCountDel + 1
							if (repCountDel == 5){
								stop(mbCountMethods[i]," unable to find ",nc," clusters when removing column ",del, 
								". \n  Stability measures and/or BSI cannot be calculated")									
							}
      					} 
        				names(clusterDel) <- rownames(objDel) 

        			if("stability"%in%validation) {
          				stabmeas <- stability(obj, Dist, del, MBcluster, clusterDel)
          				MBmeasures["APN",ind] <- MBmeasures["APN",ind] + stabmeas["APN"]
          				MBmeasures["AD",ind]  <- MBmeasures["AD",ind]  + stabmeas["AD"]
          				MBmeasures["ADM",ind] <- MBmeasures["ADM",ind] + stabmeas["ADM"]
          				MBmeasures["FOM",ind] <- MBmeasures["FOM",ind] + stabmeas["FOM"]
        				}
        			if("biological"%in%validation) {        				
						tmp <- BSI(MBcluster,clusterDel,annotation=annotation,
                     		names=rownames(objDel), category=GOcategory, goTermFreq=goTermFreq,
                     		dropEvidence=dropEvidence)
         	 			MBmeasures["BSI",ind] <- MBmeasures["BSI",ind] + tmp
        			}
        
        			if (del/ncol(obj) > 0.25 & co.del==0) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 25% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
             		 		print(paste("BSI 25% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
        			else if (del/ncol(obj) > 0.50 & co.del==1) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 50% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
              				print(paste("BSI 50% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
        			else if (del/ncol(obj) > 0.75 & co.del==2) 
          			{
            			if(clVerbose & "stability"%in%validation) 
              				print(paste("Stability validation 75% finished", mbCountMethods[i], nc, "clusters"))
            			if(clVerbose & "biological"%in%validation) 
              				print(paste("BSI 75% finished", mbCountMethods[i], nc, "clusters"))            
            			co.del <- co.del+1
          			}
      			} #END OF del LOOP
      			if(clVerbose & "stability"%in%validation)
        			print(paste("Finished stability validation", mbCountMethods[i], nc, "clusters"))
      			if(clVerbose & "biological"%in%validation)
      				print(paste("Finished BSI", mbCountMethods[i], nc, "clusters"))
    		} #END of STABILITY measures
    		ind <- ind+1  #ind tracks number clusters
  		} #END OF NC LOOP
  
  		if ("stability"%in%validation) {
    		MBmeasures["APN",] <- MBmeasures["APN",]/ncol(obj)
    		MBmeasures["AD",] <-  MBmeasures["AD",]/ncol(obj)
    		MBmeasures["ADM",] <- MBmeasures["ADM",]/ncol(obj)
    		MBmeasures["FOM",] <- MBmeasures["FOM",]/ncol(obj)
  		}
  		if ("biological"%in%validation) {
    		MBmeasures["BSI",] <- MBmeasures["BSI",]/ncol(obj)
  		}
  		allClusterObj[[i]] <- MBclusterObj
  		allMeasures[,,i] <- MBmeasures
  	} # END OF MBCOUNT METHODS LOOP	
	list(clusterObj=allClusterObj, measures=allMeasures)		
}