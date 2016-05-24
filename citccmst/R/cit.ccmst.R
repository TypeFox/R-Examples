
# Rd
# description >> assign expression data sample(s) to CIT Colon Cancer Molecular Subtype(s)
# argument
# item >> data >> a data.frame of expression data with id as rownames
# item >> data.annot >> a data.frame of data probes annotations
# item >> data.colId >> name of the column in data.annot containing data probes id
# item >> data.colMap >> name of the column in data.annot containing data probes names to map
# item >> citccmst.annot >> affymetrix annotation data.frame, if NULL (default) take the embedded annotation in object citccmst$data.annot
# item >> citccmst.colId >> name of the column in citccmst.annot corresponding to rownames of citccmst$data. Default "Probe.Set.ID"
# item >> citccmst.colMap >> name of the column in citccmst.annot containing the same annotation as in data.colMap
# item >> dist.method >> metric to compute distance to assign a sample to a subtype ("pearson", "dlda", "dqda","euclidian"). Default "dqda". 
# item >> dist.difftopcentcutoff >> cut-off on the differences between distances to centroids. If the distance is inferior to this cut-off for n centroids the sample is assigned to the n subtypes in the output variable citccmst.mixed. If NULL, the cut-off is defined as the 1st decile of the difference between the top 2 closest centroids on data used to compute centroids.
# item >> dist.disttocentcutoff >> cut-off on the mad (median absolute deviation) of distances to the centroid to define a sample as outlier. If the distance to the centroid of the assigned subtype is superior to \code{sdisttocent}*mad(distances of centroids samples to this centroid)
# item >> dist.maxcutoff >> samples for which nearest centroid is above this threshold are discarded (used only if \code{dis.meth} = "pearson" or "spearman")
# item >> plot >> if TRUE plot an acp of cit data used to classify, and of the input data with subtype affectation and dist to centroid class
# value >> a data.frame with 4 columns : "citccmst" classification to the closest of the 6 subtype centroids, "citccmst.mixed" classification to the n closest subtypes depending on dist.difftopcentcutoff, "citccmst.core" classification without outlier and mixed samples and citccmst.confidence a confidence assignation annotation (CORE, MIXED or OUTLIER)
# author >> Laetitia Marisa
# keyword >> methods
# examples >> load(list.files(system.file("extdata", package="citccmst"), full.names=TRUE))
# examples >> citvalid.exp.annot <- data.frame(id=rownames(citvalid.exp.norm), stringsAsFactors=FALSE, row.names=rownames(citvalid.exp.norm) )
# examples >> citccmst <- cit.assignCcmst(     data=citvalid.exp.norm,
# examples >>                                  data.annot=citvalid.exp.annot,
# examples >>                                  data.colId="id",
# examples >>                                  data.colMap="id" ,
# examples >>                                  citccmst.colMap="Probe.Set.ID",
# examples >>                                  plot=TRUE
# examples >>                            )
# examples >> head(citccmst)
# end

cit.assignCcmst <- function(  data,
                              data.annot,
                              data.colId="Probe.Set.ID",
                              data.colMap=c("Probe.Set.ID","Gene.Symbol","Ensembl","UniGene.ID")[1],
                              citccmst.annot=NULL,
                              citccmst.colId="Probe.Set.ID",
                              citccmst.colMap=c("Probe.Set.ID","Gene.Symbol","Ensembl","UniGene.ID")[1],
                              dist.method="dqda",
                              dist.difftopcentcutoff=NULL,
                              dist.disttocentcutoff=NULL,
                              dist.maxcutoff=NULL,
                              plot=FALSE
                              )
{

  #data("citccmst",  envir=sys.frame(sys.nframe()) ) #get citccmst

  if( is.null(citccmst.annot) ){
   citccmst.annot <- citccmst$data.annot
  }
  if( is.null(rownames(data.annot)) ){
    rownames(data.annot) <- data.annot[,data.colId]
  }
  # check that rownames == data.colId
  if( ! all(rownames(data)%in% rownames(data.annot)) )
    stop("data should have the same rownames as data.annot.")

  if(!is.data.frame(data)){
    if( is.matrix(data)){
      data <- as.data.frame(data)
    }else{
      stop("data should be either a matrix or a dataframe.")
    }
  }

    ##
    ## get cit probe annotation to map
    ##
    idcitccmst <- rownames(citccmst$data)
    mapcitccmst <- citccmst.annot[idcitccmst,citccmst.colMap]

    ## rm not pertinent probes
    torm <- which(mapcitccmst %in% c("---","", NA))
    if(length(torm)){
      idcitccmst <- idcitccmst[-torm]
      mapcitccmst <- mapcitccmst[-torm]
    }


    ##
    ## get input annotation matching cit
    ##
    iddata <- rownames(data)
    mapdata <- data.annot[iddata,data.colMap]

    map <- mapcitccmst[mapcitccmst %in% mapdata]
    names(map) <- idcitccmst[mapcitccmst %in% mapdata]
    if(length(map)==0)
      stop("None of the probes are mapped to cit data. Check column names entered.")


    map2 <- mapdata[mapdata %in% mapcitccmst]
    names(map2) <- iddata[mapdata %in% mapcitccmst]

    nmappedprobes <- table(unique(mapcitccmst) %in% unique(mapdata))["TRUE"]
    cat("Mapping - ", nmappedprobes["TRUE"], "/", length(unique(mapcitccmst))," original probes.\n", sep="")


    ## data aggregation
    citccmst.red <- cit.dfAggregate( citccmst$data[names(map),], map,  MARGIN = 1, fAggreg = function(x) median(x, na.rm=TRUE))

    data.red <- cit.dfAggregate(data[names(map2),], map2,  MARGIN = 1, fAggreg = function(x) median(x, na.rm=TRUE) )
    data.red <- data.red[rownames(citccmst.red),]


    ##
    ## Classify data
    ##

    median.na <- function(m) median(m, na.rm=T)

    grp <- citccmst$data.cl
    cc <- cit.centroids(citccmst.red, classes = grp, rowCentering = median.na)
    res <- cit.distToCentroids(centroids=cc$centroids, d=data.red, dist.meth = dist.method, maxDist = dist.maxcutoff, sdifftop=dist.difftopcentcutoff, sdisttocent=dist.disttocentcutoff)

    wccit <- round(sum(diag(table(cc$distToCentroids$predCore, grp)))/sum(table(grp))*100, 1)
  
 
    grpcit <- data.frame(row.names=names(data))
    grpcit[, "citccmst"] <- res$pred[names(data)]
    grpcit[, "citccmst.mixed"] <- res$predwmixed[names(data)]
    grpcit[, "citccmst.core"] <- res$predCore[names(data)]
    grpcit[, "citccmst.confidence"] <- res$group.confidence[names(data)]
    attributes(grpcit)$distmethod <- dist.method
    attributes(grpcit)$maxdist <- res$cutoffmaxdist
    attributes(grpcit)$nb.mapped.probes <- paste(nmappedprobes,"/",nrow(citccmst$data),sep="")
    attributes(grpcit)$citmc <-  paste(wccit, "%", sep="")
    attributes(grpcit)$scoreGroup <- c(res$cutoffdiffdist, res$cutoffdisttocent)
  
    ## plot results to check
    if(plot){
      ccmstcol <- function(x) c(C1 = "red", C2 = "blue", C3 = "green", C4 = "orange", C5 = "purple", C6 = "gray")[x]

      citpca <- prcomp(t(citccmst.red), center=TRUE, scale.=FALSE)
      
      # check for na's
      data.red.wona <- data.red[apply(is.na(data.red),1,sum)==0,]
      datapca <- prcomp(t(data.red.wona), center=TRUE, scale.=FALSE)
  
      newgrp <- grpcit$citccmst
  
       
      par(mfrow=c(2,1), mar=c(3,5,2,1))
      plot(citpca$x[,c(1:2)], col = ccmstcol(grp), pch=16, asp=1, main="PCA on CIT data used to classify dataset")
      for (g in unique(grp)){
        wg <- which(grp == g)
        if(length(wg)>1){
        mu <- c(mean(citpca$x[wg, 1]), mean(citpca$x[wg, 2]))
        for(j in wg)
          lines(rbind(mu, as.numeric(citpca$x[j,1:2])), col = ccmstcol(g), lty = "dotted")
        }
        legend("bottomright", legend = levels(as.factor(grp)), fill = ccmstcol(sort(unique(grp))), bty = "n", cex = 1)
      }
  

  
      vpch <- rep(16, ncol(data))
      vpch[is.na(grpcit$citccmst.core)] <- 1
      plot(datapca$x[,c(1:2)], col = ccmstcol(newgrp), pch=vpch, asp=1, main="Input dataset")
      for (g in unique(newgrp)){
        wg <- which(newgrp == g)
        if(length(wg)>1){
        mu <- c(mean(datapca$x[wg, 1]), mean(datapca$x[wg, 2]))
        for(j in wg)
          lines(rbind(mu, as.numeric(datapca$x[j,1:2])), col = ccmstcol(g), lty = "dotted")
        }
        legend("bottomright", legend = levels(as.factor(newgrp)), fill = ccmstcol(sort(unique(newgrp))), bty = "n", cex = 0.8)
        legend("topright", legend = c( "core","uncertain"), pch=c(16,1), bty = "n", cex = 0.8, title="sample classification")
      }
  
  
    }
  
    grpcit
  
}



#cit.oncotypeLikeClassifier <- function(){
#                              
#                              
#}


