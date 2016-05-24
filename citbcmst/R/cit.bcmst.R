
# Rd
# description >> assign expression data sample(s) to CIT Breast Cancer Molecular Subtype(s)
# argument
# item >> data >> a data.frame of expression data with id as rownames
# item >> data.annot >> a data.frame of data probes annotations
# item >> data.colId >> name of the column in data.annot containing data probes id
# item >> data.colMap >> name of the column in data.annot containing data probes names to map
# item >> citbcmst.annot >> affymetrix annotation data.frame, if NULL (default) take the embeded annotation in object citbcmst$data.annot
# item >> citbcmst.colId >> name of the column in citbcmst.annot corresponding to rownames of citbcmst$data. Default "Probe.Set.ID"
# item >> citbcmst.colMap >> name of the column in citbcmst.annot containing the same annotation as in data.colMap
# item >> dist.method >> metric to compute distance to assign a sample to a subtype. Default for Affymetrix data "dlda". For other platforms "pearson".
# item >> dist.difftopcentcutoff >> cut-off on the differences between distances to centroids. If the distance is inferior to this cut-off for n centroids the sample is assigned to the n subtypes in the output variable citbcmst.mixed. If NULL, the cut-off is defined as the 1st decile of the difference between the top 2 closest centroids on data used to compute centroids.
# item >> dist.disttocentcutoff >> cut-off on the mad (median absolute deviation) of distances to the centroid to define a sample as outlier. If the distance to the centroid of the assigned subtype is superior to \code{sdisttocent}*mad(distances of centroids samples to this centroid)
# item >> dist.maxcutoff >> samples for which nearest centroid is above this threshold are discarted (used only if \code{dis.meth} = "pearson" or "spearman")
# item >> plot >> if TRUE plot an acp of citdata used to classify, and of the input data with subtype affectation and dist to centroid class
# value >> a data.frame with 4 columns : "citbcmst" classification to the closest of the 6 subtypes, "citbcmst.mixed" classification to the n closest subtypes depending on dist.difftopcentcutoff, "citbcmst.core" classification without outlier and mixed samples and citbcmst.confidence a confidence assignation annotation (CORE, MIXED or OUTLIER)
# author >> Laetitia Marisa
# keyword >> methods
# examples >> load(list.files(system.file("extdata", package="citbcmst"), full.names=TRUE)[1])# load exp.norm.bertheau07 object stored in /inst/extdata
# examples >> exp.annot.bertheau07 <- data.frame(id=rownames(exp.norm.bertheau07), stringsAsFactors=FALSE, row.names=rownames(exp.norm.bertheau07) )
# examples >> citbcmst.bertheau07 <- cit.assignBcmst(   data=exp.norm.bertheau07,
# examples >>                                           data.annot=exp.annot.bertheau07,
# examples >>                                           data.colId="id",
# examples >>                                           data.colMap="id" ,
# examples >>                                           citbcmst.colMap="Probe.Set.ID",
# examples >>                                           dist.method="dlda",
# examples >>                                           plot=TRUE
# examples >>                                       )
# end

cit.assignBcmst <- function(  data,
                              data.annot,
                              data.colId="Probe.Set.ID",
                              data.colMap=c("Probe.Set.ID","Gene.Symbol","Ensembl","UniGene.ID")[1],
                              citbcmst.annot=NULL,
                              citbcmst.colId="Probe.Set.ID",
                              citbcmst.colMap=c("Probe.Set.ID","Gene.Symbol","Ensembl","UniGene.ID")[1],
                              dist.method="dlda",
                              dist.difftopcentcutoff=NULL,
                              dist.disttocentcutoff=NULL,
                              dist.maxcutoff=NULL,
                              plot=FALSE
                              )
{

  #warning("This classifier is optimised for Affymetrix (U133plus2) data.")

  data("citbcmst",  envir=sys.frame(sys.nframe()) ) #get citbcmst

  if( is.null(citbcmst.annot) ){
   citbcmst.annot <- citbcmst$data.annot
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
    idcitbcmst <- rownames(citbcmst$data)
    mapcitbcmst <- citbcmst.annot[idcitbcmst,citbcmst.colMap]

    ## rm not pertinent probes
    torm <- which(mapcitbcmst %in% c("---","", NA))
    if(length(torm)){
      idcitbcmst <- idcitbcmst[-torm]
      mapcitbcmst <- mapcitbcmst[-torm]
    }


    ##
    ## get input annotation matching cit
    ##
    iddata <- rownames(data)
    mapdata <- data.annot[iddata,data.colMap]

    map <- mapcitbcmst[mapcitbcmst %in% mapdata]
    names(map) <- idcitbcmst[mapcitbcmst %in% mapdata]
    if(length(map)==0)
      stop("None of the probes are mapped to cit data. Check column names entered.")


    map2 <- mapdata[mapdata %in% mapcitbcmst]
    names(map2) <- iddata[mapdata %in% mapcitbcmst]

#    nmappedprobes <- table(mapcitbcmst %in% mapdata)["TRUE"]
#    cat("Mapping - ", nmappedprobes["TRUE"], "/", nrow(citbcmst$data)," original probes.\n", sep="")
    nmappedprobes <- table(unique(mapcitbcmst) %in% unique(mapdata))["TRUE"]
    cat("Mapping - ", nmappedprobes["TRUE"], "/", length(unique(mapcitbcmst))," original probes.\n", sep="")


    ## data aggregation
    citbcmst.red <- cit.dfAggregate( citbcmst$data[names(map),], map,  MARGIN = 1, fAggreg = function(x) median(x, na.rm=TRUE))

    data.red <- cit.dfAggregate(data[names(map2),], map2,  MARGIN = 1, fAggreg = function(x) median(x, na.rm=TRUE) )
    data.red <- data.red[rownames(citbcmst.red),]


    ##
    ## Classify data
    ##
    grp <- c(citbcmst$data.cl,rep(NA,ncol(data.red)))
    res <- cit.centroids(cbind(citbcmst.red, data.red), classes = grp, rowCentering = NA, dist.meth = dist.method, maxDist = dist.maxcutoff, sdifftop=dist.difftopcentcutoff, sdisttocent=dist.disttocentcutoff)
  
    #wccit <- round(sum(diag(table(res$distToCentroids$pred, grp)))/sum(table(grp))*100, 1)
    wccit <- round(sum(diag(table(res$distToCentroids$predCore, grp)))/sum(table(grp))*100, 1)
  
#    if(getOption("verbose")==TRUE){
#      table(res$distToCentroids$pred,c(rep("CIT",ncol(citbcmst.red)),rep("input",ncol(data.red))))[c("basL","mApo","lumC","lumB","lumA","normL"),]
#    }
    cat("Classification - ",wccit,"% of cit data well classified after data reduction.\n", sep="")
  
    cit.bcmstGroups <- function( bcmst )
    {
      bcmstgrp <- c(
      "basL" = "basL",
      "mApo" = "mApo",
      "lumA" = "lumA",
      "lumB" = "lumB",
      "lumC" = "lumC",
      "normL" = "normL",
      "basLmApo"="basLmApo",
      "basLlumC"="basLlumC",
      "basLlumB"="basLlumB",
      "basLlumA"="basLlumA",
      "basLnormL"="nbasL",
      "mApobasL"="basLmApo",
      "mApolumC"="mApolumC",
      "mApolumB"="mApolumB",
      "mApolumA"="mApolumA",
      "mAponormL"="nmApo",
      "lumAlumB"="lumAB",
      "lumAlumC"="lumAC",
      "lumAnormL"="nlumA",
      "lumAbasL"="basLlumA",
      "lumAmApo"="mApolumA",
      "lumBlumA"="lumAB",
      "lumBlumC"="lumBC",
      "lumBnormL"="nlumB",
      "lumBbasL"="basLlumB",
      "lumBmApo"="mApolumB",
      "lumClumA"="lumAC",
      "lumClumB"="lumBC",
      "lumCnormL"="nlumC",
      "lumCbasL"="basLlumC",
      "lumCmApo"="mApolumC",
      "normLlumA"="nlumA",
      "normLlumB"="nlumB",
      "normLlumC"="nlumC",
      "normLbasL"="nbasL",
      "normLmApo"="nmApo",
      "lumAlumBlumC"="lumABC",
      "lumAlumBnormL"="nlumAB",
      "lumAlumCnormL"="nlumAC",
      "lumBlumCnormL"="nlumBC",
      "lumCmAponormL"="nmApolumC",
      "lumAlumBlumCnormL"="nlumABC",
      "lumBlumCmApo"="mApolumBC"
    )
    v <- as.vector(bcmstgrp[as.character(bcmst)])
    v[is.na(v)] <- as.character(bcmst)[is.na(v)]
    v
    }
  
    grpcit <- data.frame(row.names=names(data))
    grpcit[, "citbcmst"] <- res$distToCentroids$pred[names(data)]
    grpcit[, "citbcmst.mixed"] <- cit.bcmstGroups(res$distToCentroids$predwmixed[names(data)])
    grpcit[, "citbcmst.core"] <- res$distToCentroids$predCore[names(data)]
    grpcit[, "citbcmst.confidence"] <- res$distToCentroids$group.confidence[names(data)]
    attributes(grpcit)$distmethod <- dist.method
    attributes(grpcit)$maxdist <- res$cutoffmaxdist
    attributes(grpcit)$nb.mapped.probes <- paste(nmappedprobes,"/",nrow(citbcmst$data),sep="")
    attributes(grpcit)$citmc <-  paste(wccit, "%", sep="")
    attributes(grpcit)$scoreGroup <- c(res$distToCentroids$cutoffdiffdist, res$distToCentroids$cutoffdisttocent)
  
    ## plot results to check
    if(plot){
      bcmstcol <- function(x) c(lumBC = "maroon3", nlumA = "lightseagreen", basL = "red", mApo = "orange", lumA = "blue", lumB = "lightskyblue", lumC = "hotpink", normL = "green3")[x]
  
      grpb <- grp[-which(is.na(grp))]
      citpca <- prcomp(t(citbcmst.red[,names(grpb)]), center=TRUE, scale.=FALSE)
      
      # check for na's
      data.red.wona <- data.red[apply(is.na(data.red),1,sum)==0,]
      datapca <- prcomp(t(data.red.wona), center=TRUE, scale.=FALSE)
  
      newgrp <- grpcit$citbcmst
  
      if( any(newgrp=="basL") ){
        signb <- sign(apply(citpca$x[grpb=="basL",c(1:2)],2,mean,na.rm=TRUE))
        signb2 <- sign(apply(datapca$x[,c(1:2)][newgrp=="basL",],2,mean,na.rm=TRUE))
        if(signb2[1] != signb[1]){
          datapca$x[,1] <- (-1)*datapca$x[,1]
        }
        if(signb2[2] != signb[2]){
          datapca$x[,2] <- (-1)*datapca$x[,2]
        }
      }else{
      if( any(newgrp=="mApo") ){
        signm <- sign(apply(citpca$x[grpb=="mApo",c(1:2)],2,mean,na.rm=TRUE))
        signm2 <- sign(apply(datapca$x[newgrp=="mApo",c(1:2)],2,mean,na.rm=TRUE))
        if(signm2[1] != signm[1]){
          datapca$x[,1] <- (-1)*datapca$x[,1]
        }
        if(signm2[2] != signm[2]){
          datapca$x[,2] <- (-1)*datapca$x[,2]
        }
      }
      }
      
      par(mfrow=c(2,1), mar=c(3,5,2,1))
      plot(citpca$x[,c(1:2)], col = bcmstcol(grpb), pch=16, asp=1, main="PCA on CIT data used to classify dataset")
      for (g in unique(grpb)){
        wg <- which(grpb == g)
        if(length(wg)>1){
        mu <- c(mean(citpca$x[wg, 1]), mean(citpca$x[wg, 2]))
        for(j in wg)
          lines(rbind(mu, as.numeric(citpca$x[j,1:2])), col = bcmstcol(g), lty = "dotted")
        }
        legend("bottomright", legend = levels(as.factor(grpb)), fill = bcmstcol(sort(unique(grpb))), bty = "n", cex = 1)
      }
  

  
      vpch <- rep(16, ncol(data))
      vpch[is.na(grpcit$citbcmst.core)] <- 1
      plot(datapca$x[,c(1:2)], col = bcmstcol(newgrp), pch=vpch, asp=1, main="Input dataset")
      for (g in unique(newgrp)){
        wg <- which(newgrp == g)
        if(length(wg)>1){
        mu <- c(mean(datapca$x[wg, 1]), mean(datapca$x[wg, 2]))
        for(j in wg)
          lines(rbind(mu, as.numeric(datapca$x[j,1:2])), col = bcmstcol(g), lty = "dotted")
        }
        legend("bottomright", legend = levels(as.factor(newgrp)), fill = bcmstcol(sort(unique(newgrp))), bty = "n", cex = 0.8)
        legend("topright", legend = c( "core","uncertain"), pch=c(16,1), bty = "n", cex = 0.8, title="sample classification")
      }
  
  
    }
  
    grpcit
  
}






