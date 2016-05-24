# Rd
# description >> define centroids for a given partition of individuals (columns) and calculate distance between individuals and centroids
# argument
# item >> d >> a \code{data.frame} of numeric data
# item >> classes >> a vector defining a partition of \code{data} columns (\code{NA} values accepted)
# item >> rowCentering >>  \code{NA}: no row centering; otherwise: function to be used for row centering
# item >> rowClassesForAggregation >>  partition of the rows , used to aggregate rows from the same class (1 aggregated row per class will be calculated)
# item >> rowClassesToKeep >> to restrict the centroid's calculation to the some aggregated rows  (= row classes)
# item >> dist.meth >> distance method used to calculate distance between individuals (columns in \code{d}) and centroids
# item >> maxDist >> individuals for which nearest centroid is above this threshold are discarted (used only if \code{dis.meth} = \code{"pearson"} or \code{"spearman"})
# item >> ... >> parameters from cit.distToCentroids function
# value >> a \code{list} with two objects : \code{centroids} and \code{dist2Centroids}
# author >> Aurelien de Reynies
# keyword >> internal
# end
cit.centroids <- function(  d,
                            classes,
                            rowCentering=c(NA,function(x)mean(x, na.rm=TRUE),function(x)median(x, na.rm=TRUE))[[3]],
                            rowClassesForAggregation=NULL,
                            rowClassesToKeep=NULL,
                            dist.meth=c("spearman","euclidian","maximum","manhattan","canberra","binary","minkowski","pearson","dlda","dqda"),
                            maxDist=0.5,
                            ...
                            )
{

    n <- length(classes)
    N <- dim(d)[2]
    if(n != N)  stop("Error - function cit.centroids: size of classes doesn't correspond to number of columns in data")
    dist.meth <- match.arg(dist.meth)

    if(!is.null(rowClassesForAggregation)) {
          d <- cit.dfAggregate(d,rowClassesForAggregation)
          if(!is.null(rowClassesToKeep)) {
               rowClassesToKeep <- intersect(rowClassesToKeep,rownames(d))
               if(length(rowClassesToKeep) < 2) stop("Error - function cit.centroids  : less than 2 rowClassesToKeep are in d")
               d <- d[rowClassesToKeep,]
          }
    }

    if(is.function(rowCentering))  d <- sweep(d,1,apply(d,1,rowCentering))

    w <- which(!is.na(classes))

    if(length(w)==0) stop("Error - function cit.centroids: classes are undefined (NA)")


    #classes <- classes[w]
    nc <- table(classes)

    m <- cit.dfAggregate(d[,w],classes[w],MARGIN=2,fAggreg=function(x) mean(x, na.rm=TRUE) )
    v <- cit.dfAggregate(d[,w],classes[w],MARGIN=2,fAggreg=function(x) var(x, na.rm=TRUE))
    va <-  apply(v,1,function(vg)  sum((nc-1)*vg)/(length(w) - length(nc)))
    vg <- apply(d[,w],1,function(x) var(x, na.rm=TRUE))

    L <- list("mean"=m,
              "var"=v,
              "aggregated var"=va,
              "global var"=vg,
              "rowCentering"=rowCentering,
              "rowClassesForAggregation"=rowClassesForAggregation,
              "rowClassesToKeep"=rowClassesToKeep,
              "centroidsdata"=list("samples"=data.frame("samplename"=colnames(d[,w]), "class"=classes[w], stringsAsFactors=FALSE),"data"=d[,w])
              )

    list("centroids"=L,
         "distToCentroids"=cit.distToCentroids(d,L,dist.meth=dist.meth,maxDist=maxDist,d.isPretreated=TRUE, ...))

}

# Rd
# description >> calculate distance between individuals and centroids; used in \code{cit.centroids()}
# argument
# item >> d >> a \code{data.frame} of numeric data
# item >> centroids >> object obtained via a call to \code{cit.centroids()} (see slot 'centroids')
# item >> dist.meth >> distance method used to calculate distance between individuals (columns in \code{d}) and centroids
# item >> maxDist >> individuals for which nearest centroid is above this threshold are discarted (used only if \code{dis.meth} = \code{"pearson"} or \code{"spearman"})
# item >> d.isPretreated >> indicate wether \code{d} is pretreated (row aggregation, row centering)...for internal use
# item >> sdifftop >> cut-off on the diffences between distances to centroids. If the distance is inferior to this cut-off for n centroids the sample is assigned to the n groups in the output variable predmixed. If NULL, the cut-off is defined as the 1st decile of the difference between the top 2 closest centroids on data used to compute centroids.
# item >> sdisttocent >> cut-off on the mad (median absolute deviation) of distances to the centroid to define a sample as outlier. If the distance to the centroid of the assigned group is superior to \code{sdisttocent}*mad(distances of centroids samples to this centroid)
# item >> verbose >> boolean if output should be displayed
# value >> a \code{list} giving the distance between individuals and centroids and the predicted class for each individuals
# author >> Aurelien de Reynies, Laetitia Marisa
# keyword >> internal
# end
cit.distToCentroids <- function(d,
                                centroids,
                                dist.meth=c("spearman","euclidian","maximum","manhattan","canberra","binary","minkowski","pearson","dlda","dqda"),
                                maxDist=0.5,
                                d.isPretreated=FALSE,
                                sdifftop=NULL,
                                sdisttocent=NULL,
                                verbose=FALSE
                                )
{

    dist.meth <- match.arg(dist.meth)

    srcs <- colnames(d)

    ## Compute distance to centroids
    if(!d.isPretreated){
          if(!is.null(centroids$rowClassesForAggregation)) {
                d <- cit.dfAggregate(d,centroids$rowClassesForAggregation)
                if(!is.null(centroids$rowClassesToKeep)) {
                     centroids$rowClassesToKeep <- intersect(centroids$rowClassesToKeep,rownames(d))
                     if(length(centroids$rowClassesToKeep) < 2) stop("Error - function cit.distToCentroids: less than 2 rowClassesToKeep are in d")
                     d <- d[centroids$rowClassesToKeep,]
                }
          }

          if(is.function(centroids$rowCentering))  d <- sweep(d,1,apply(d,1,centroids$rowCentering))
    }

    ## add centroid data to compute distance on centroids data
    dc <- NULL
    if( !is.null(centroids$centroidsdata) ){  # for previous version of cit.centroids
      if(!all(centroids$centroidsdata$samples$samplename %in% colnames(d)) ){
          dc <- centroids$centroidsdata$data
          colnames(dc) <- paste("centroid.",colnames(dc),sep="")
          d <- cbind( d, dc )
      }else{
          colnames(d)[ colnames(d) %in% centroids$centroidsdata$samples$samplename] <-  paste("centroid.",colnames(d)[ colnames(d) %in% centroids$centroidsdata$samples$samplename],sep="")
      }
    }


    N <- ncol(d)
    n <- ncol(centroids$"mean")

    if(dist.meth %in% c("dlda","dqda")){
       sumlogvar <- apply(log(centroids$"var"),2,sum, na.rm=TRUE)

       if(dist.meth == "dlda") {
            tmp <- apply(d,2,function(z) apply((z-centroids$"mean")^2/centroids$"aggregated var",2,sum, na.rm=TRUE))
            tdist <- as.data.frame(t(tmp))
       }
       if(dist.meth == "dqda"){
            tmp <-  apply(d,2,function(z) apply((z-centroids$"mean")^2/centroids$"var",2,sum, na.rm=TRUE)+sumlogvar)
            tdist <- as.data.frame(t(tmp))
       }
    }else{
       d2 <- t(cbind(d,centroids$"mean"))

       tdist <- as.matrix(cit.dist(d2 ,meth= dist.meth,diag=TRUE))
       tdist <- as.data.frame(tdist[1:N,(N+1):(N+n)])
    }


    rownames(tdist) <- names(d)
    names(tdist) <- names(centroids$"mean")

    ## affect group to the closest centroids
    pred  <- apply(tdist,1,function(z) names(centroids$"mean")[which.min(z)])
    mind  <- apply(tdist,1,min)


    ## define mixed samples, i.e. samples between n groups

    difftofirst <- function(x){
      m <- min(x)
      x-m
    }

    difftop <- apply( tdist, 1, function(x) { diff(sort(x)[1:2]) })

    if(is.null(sdifftop) ){
      if( length(grep("centroid.",names(difftop),value=TRUE))>0)
        sdifftop <- round(quantile(difftop[grep("centroid.",names(difftop),value=TRUE)], 0.01,na.rm=TRUE),3)
      else
        sdifftop <- round(quantile(difftop, 0.01,na.rm=TRUE),3)
      if(verbose==TRUE)
        cat("Estimated cut-off for difference between distance to centroids = ", sdifftop, "\n", sep="")
    }
    predf <- sapply(1:nrow(tdist), function(i){ if(difftop[i]<=sdifftop){
                                                            w <- which(difftofirst(tdist[i,])<=sdifftop)
                                                            paste(sort(colnames(tdist[i,])[w]), collapse="")
                                                          }else{
                                                            colnames(tdist[i,])[which.min(tdist[i,])]
                                                          }})
    names(predf) <- names(pred)


    ## define outlier samples, i.e. samples outside centroid limits
    ## ------------------------------------------------------------

    if( length(grep("centroid.",names(difftop),value=TRUE))>0 ){

      sam <- grep("centroid.",names(difftop),value=TRUE)

      wsure <- sam[ centroids$centroidsdata$samples[match(sub("centroid.","",sam), centroids$centroidsdata$samples$samplename),"class"] == predf[sam] ]

      if(length(wsure)==0)
        stop("No concordance between pred and predf on centroids data. sdifftop is too high.\n")

      coresettab <- data.frame(row.names=wsure)
      coresettab$groups <-  predf[wsure]
      coresettab <- cbind( coresettab, as.matrix(tdist[wsure,]))
      coresettab$disttocent <-   mind[wsure]
      coresettab$difftop <- apply( tdist[wsure,], 1, function(x) { diff(sort(x)[1:2]) })

      refcoreset <- NULL
      for( g in names(tdist) ){
        L <- split(coresettab[,g], coresettab[,"groups"]==g)["TRUE"]
        inf <- lapply(L, function(x) c(median(x),max(x),mad(x)) )
        refcoreset <- cbind( refcoreset, matrix(unlist(inf), ncol=1, dimnames=list(c("med","max","mad"),g)))
      }

      if(is.null(sdisttocent)){
        sdisttocent <- max(round((refcoreset[2,]-refcoreset[1,])/refcoreset[3,]))
        if(verbose==TRUE)
          cat("Estimated cut-off on distance to centroids mad  =", sdisttocent, "\n")
        sdisttocent <- refcoreset["med",pred]+sdisttocent*refcoreset["mad",pred]
      }else{
        if(length(sdisttocent)==1)
           sdisttocent <- refcoreset["med",pred]+sdisttocent*refcoreset["mad",pred]
      }

    if(verbose){
      print(refcoreset)
      print(unique(cbind(pred, round(sdisttocent,3))))
    }

    scoregroup <- c("TRUE"="OUTLIER","FALSE"="CORE")[as.character(mind>sdisttocent)]
    scoregroup[which(difftop <= sdifftop)] <- "MIXED"
    names(scoregroup) <- names(difftop)

    pred2 <- pred
    pred2[ ! scoregroup %in% "CORE"] <- NA

    if( !is.null(maxDist) & dist.meth %in% c("pearson","spearman")){
      pred2[which(mind > maxDist)] <- NA
      scoregroup[which(mind > maxDist)] <- "OUTLIER"
    }else{
      maxDist <- NULL
    }

    }else{ # if no centroid data
      scoregroup <- rep("ND", length(pred))
      scoregroup[which(difftop <= sdifftop)] <- "MIXED"
      pred2 <- pred
      pred2[ scoregroup %in%"MIXED"] <- NA

      if( !is.null(maxDist) & dist.meth %in% c("pearson","spearman")){
         pred2[which(mind > maxDist)] <- NA
         scoregroup[which(mind > maxDist)] <- "OUTLIER"
      }else{
         pred2[which(mind > quantile(mind,0.95))] <- NA
         scoregroup[which(mind > quantile(mind,0.95))] <- "OUTLIER"
         maxDist <- quantile(mind,0.95)
      }
    }


    if( is.null(dc) ){
      rownames(tdist) <- sub("centroid.","",rownames(tdist))
      names(difftop) <- sub("centroid.","",names(difftop))
      names(mind) <- sub("centroid.","",names(mind))
      names(pred) <- sub("centroid.","",names(pred))
      names(predf) <- sub("centroid.","",names(predf))
      names(pred2) <- sub("centroid.","",names(pred2))
      names(scoregroup) <- sub("centroid.","",names(scoregroup))
    }else{
      tdist <- tdist[srcs,]
      difftop <- difftop[srcs]
      mind <- mind[srcs]
      pred <- pred[srcs]
      predf <- predf[srcs]
      pred2 <- pred2[srcs]
      scoregroup <- scoregroup[srcs]
    }



    list("dist.scores"            = tdist,
         "distToNearestCentroid"  = mind,
         "diffDistTopCentroids"   = difftop,
         "pred"                   = pred,
         "predwmixed"             = predf,
         "predCore"               = pred2,
         "group.confidence"       = scoregroup,
         "dist.meth"              = dist.meth,
         "cutoffdiffdist"         = sdisttocent,
         "cutoffdisttocent"       = sdifftop,
         "cutoffdistmax"          = maxDist
         )

}




# Rd
# description >> aggregate each column of \code{d} by \code{byvar} using function \code{FUN}
# argument
# item >> data >> a dataframe of numeric data
# item >> partition >> a vector on which to aggregate
# item >> MARGIN >>  1 : to aggregate rows ; 2 : to aggregate columns
# item >> fAggreg >> the function to aggregate \code{data} by \code{partition}
# value >> the reduced data.frame \code{data}
# author >> Aurelien de Reynies
# keyword >> internal
# end
cit.dfAggregate <- function( data, partition, MARGIN=1, fAggreg=function(x){mean(x, na.rm=TRUE)}){

    cMARGIN <- setdiff(c(1,2),MARGIN)
    n <- length(partition)
    N <- dim(data)[MARGIN]

    p <- dim(data)[cMARGIN]

    if(n != N)  stop("Error - function cit.dfAggregate : size of partition doesn't correspond to data dimension")
    l <- split(1:N,partition)

    d <- data
    if(MARGIN == 2) d <- t(data)

    d <- matrix(sapply(l,function(i)
                                if(length(i)==1){unlist(d[i,])}
                                else{ apply(d[i,],2,fAggreg)}
                ), ncol=p, byrow=TRUE)

    d <- as.data.frame(d)
    rownames(d) <- names(l)
    names(d)    <- dimnames(data)[[cMARGIN]]

    if(MARGIN == 2) d <- as.data.frame(t(d))

    d
}


# Rd
# description >> This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
# argument
# item >> x >> a numeric matrix, data frame or "dist" object.
# item >> meth >> the distance measure to be used. This must be one of "pearson","spearman","euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
# item >> use >>  method for computing covariance when NAs, choice between "pairwise.complete.obs", "all.obs" and "complete.obs" (cf cor function parameter)
# item >> diag >> logical value indicating whether the diagonal of the distance matrix should be printed
# item >> upper >> logical value indicating whether the upper triangle of the distance matrix should be printed
# item >> p >> the power of the Minkowski distance.
# item >> replaceNA >> a boolean indicating if NA value should be replace by
# value >> required
# author >> Aurelien de Reynies, Mickael Guedj
# keyword >> internal
# seealso >> \link[stats]{dist}
# end
cit.dist <- function(x,
                  meth = "pearson",
                  use = "pairwise.complete.obs",
                  diag = FALSE,
                  upper = FALSE,
                  p = 2,
                  replaceNA = TRUE) {
    res <- NULL

    DIST.METHODS <- c("euclidian",
                   "maximum",
                   "manhattan",
                   "canberra",
                   "binary",
                   "minkowski",
                   "pearson",
                   "spearman")

    is.ok <- function(obj,fun=all){
      if(is.null(obj)) return(FALSE)
      if(length(obj)==0)return(FALSE)
      fun(!is.na(obj) & !is.nan(obj) & !is.infinite(obj))
    }


  if(!is.na(meth) & !is.null(meth) & meth %in% DIST.METHODS){
    if(meth %in% DIST.METHODS[1:6]){
      res <- dist(x,method=meth,diag=diag,upper=upper,p=p)
      maxdist <- ceiling(max(res,na.rm=TRUE))
    }else{
      maxdist <- 2
      if(!is.ok(use))   { use <- "all.obs" ; print("NB - function cit.dist (utile clustering) : parameter 'use' was set to default value = 'all.obs'") }                           # modified by MG (01.08.08)
      res <- as.dist(1-cor(t(x),use=use,method=meth),diag=diag,upper=upper)
    }
  }
  wNA <- which(is.na(res))
  if(length(wNA)>0){
    if(replaceNA){
      res[wNA] <- maxdist
      print(paste("Warning - function cit.dist (utile clustering) : d(x,y) = NA is replaced by 0 if x=y,",maxdist,"otherwise"))
      print("NB : maybe an other distance metric (ex. euclidian, manhattan,...) could avoid inconsistencies !!")
    }else{
      print("Warning - function cit.dist (utile clustering) : there are NA values in the distance matrix")
    }
  }
  if(is.null(res)) stop("ERROR - function cit.dist : uncorrect value for parameter 'meth' and/or 'use'")
  return(res)
}





