
## Defines global variables for package checks
if(getRversion() >= "2.15.1")  globalVariables(c("biocLite"))


#####################################################################################
## clValid Functions
#####################################################################################

clValid <- function(obj, nClust, clMethods="hierarchical", validation="stability", maxitems=600,
                    metric="euclidean", method="average", neighbSize=10,
                    annotation=NULL, GOcategory="all", goTermFreq=0.05,
                    dropEvidence=NULL, verbose=FALSE, ...) {



  clMethods <- tolower(clMethods)
  clMethods <- match.arg(clMethods,c("hierarchical","kmeans","diana","fanny","som","model","sota","pam","clara","agnes"), several.ok=TRUE)

  if("som"%in%clMethods) {
    if(!require(kohonen)) {
      stop("package 'kohonen' required for clustering using SOM")
    }
  }

  if("model"%in%clMethods) {
    if(!require(mclust)) {
      stop("package 'mclust' required for model-based clustering")
    }
  }

  validation <- match.arg(validation,c("stability","internal","biological"),several.ok=TRUE)
  metric <- match.arg(metric,c("euclidean", "correlation", "manhattan")) ## used for hierarchical, diana, fanny, agnes, pam
  method <- match.arg(method,c("ward", "single", "complete", "average")) ## for hclust, agnes
  GOcategory <- match.arg(GOcategory, c("all","BP","CC","MF"))

  switch(class(obj),
         matrix = mat <- obj,
         ExpressionSet = mat <- Biobase::exprs(obj),
         data.frame = {
           if(any(!sapply(obj,class)%in%c("numeric","integer")))
             stop("data frame 'obj' contains non-numeric data")
           mat <- as.matrix(obj)
         },
         ##         dist = Dist <- obj,
         stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet object"))

  if (nrow(mat)>maxitems) {
    if (interactive()) {
      cat("\nThe number of items to be clustered is larger than 'maxitems'\n")
      cat("The memory and time required may be excessive, do you wish to continue?\n")
      cat("(y to continue, any other character to exit) ")
      ans <- tolower(substr(readLines(n=1),1,1))
      if (ans!="y") {
        stop("Exiting clValid, number of items exceeds 'maxitems'")
      }
    } else {
      stop("The number of items to be clustered is larger than 'maxitems'\n  Either decrease the number of rows (items) or increase 'maxitems'\n")
    }
  }


  if ("clara"%in%clMethods & metric=="correlation")
    warning("'clara' currently only works with 'euclidean' or 'manhattan' metrics - metric will be changed to 'euclidean'  ")

  if(is.character(annotation) & length(grep(".db", annotation))==0) {
    annotation <- paste(annotation, ".db", sep="")
  }
  ## Convert annotation list to annotation table

  if (is.list(annotation)) {
    if(is.null(rownames(mat))) {
      stop("rownames of data must be present to specify biological annotation from file")
    }
    annotation <- annotationListToMatrix(annotation, genenames=rownames(mat))
  }

  if ("biological"%in%validation & is.null(annotation)) {
    stop("annotation must be specified in order to use biological validation")
  }
  if ("biological"%in%validation & is.character(annotation)) {
    if(!require(Biobase) | !require(GO.db) | !require(annotate)) {
      stop("packages 'Biobase', 'GO.db', and 'annotate' required for 2nd type of biological validation \n
these can be downloaded from Bioconductor (www.bioconductor.org)")
    }
  }
  if("biological"%in%validation & is.null(rownames(mat))) {
    stop("rownames of data must be present to use biological validation")
  }

#  if (!is.matrix(mat) | !is.numeric(mat))
#    stop("argument 'mat' must be a numeric matrix")

  nClust <- floor(nClust)
  if (any(nClust<1))
    stop("argument 'nClust' must be a positive integer vector")

  if(metric=="correlation")
    Dist <- as.dist(1-cor(t(mat), use="pairwise.complete.obs"))  else
  Dist <- dist(mat,method=metric)

  clusterObjs <- vector("list",length(clMethods))
  names(clusterObjs) <- clMethods

  measures <- c(if("stability"%in%validation) c("APN","AD","ADM","FOM"),
                if("internal"%in%validation) c("Connectivity","Dunn","Silhouette"),
                if("biological"%in%validation) c("BHI","BSI"))
  validMeasures <- array(dim=c(length(measures),length(nClust),length(clMethods)))
  dimnames(validMeasures) <- list(measures,nClust,clMethods)

  for (i in 1:length(clMethods)) {

    cvalid <- vClusters(mat,clMethods[i],nClust, validation=validation,
                        Dist=Dist, method=method, metric=metric, annotation=annotation,
                        GOcategory=GOcategory, goTermFreq=goTermFreq, neighbSize=neighbSize,
                        dropEvidence=dropEvidence, verbose=verbose, ...)
    clusterObjs[[i]] <- cvalid$clusterObj
    validMeasures[,,i] <- cvalid$measures
  }

  if(is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
    warning("rownames for data not specified, using 1:nrow(data)")
  }

  new("clValid", clusterObjs=clusterObjs, measures=validMeasures, measNames=measures,
      clMethods=clMethods, labels=rownames(mat), nClust=nClust, validation=validation,
      metric=metric,method=method, neighbSize=neighbSize,  GOcategory=GOcategory,
      goTermFreq=goTermFreq, annotation=annotation,
      call=match.call())
}




#####################################################################
## Functions for Validation Measures
## Stability, Internal, and Biological
#####################################################################

#####################################################################
## Stability Measures
####################################################################
## APN, AD, ADM, and FOM
## All measures in [0,infty] and should be minimized
#####################################################################


stability <- function(mat, Dist=NULL, del, cluster, clusterDel, method="euclidean") {

  any.na <- any(is.na(mat))
  obsNum <- 1:nrow(mat)
  nc1 <- length(table(cluster))
  nc2 <- length(table(clusterDel))
  stabmeas <- numeric(4)
  names(stabmeas) <- c("APN","AD","ADM","FOM")

  ## measure APN
  ## calculate a ncxnc matrix of proportion of non-overlaps in the two collection of nc clusters
  overlap <- xtabs(~cluster + clusterDel)
  ## measure AD
  ## calculate a ncxnc matrix of average-distance in the two collection of nc clusters
  dij <- matrix(rep(NA,nc1*nc2),nc1,nc2)

  if (is.null(Dist)) matDist <- as.matrix(dist(mat, method=method))
  if (class(Dist)=="dist") matDist <- as.matrix(Dist)
  if (class(Dist)=="matrix") matDist <- Dist

  ## measure ADM
  ## calculate a ncxnc matrix of distance-average in the two collection of nc clusters
  dij2 <- matrix(rep(NA,nc1*nc2),nc1,nc2)
  ii <- 1
  for (i in sort(unique(cluster))) {
    jj <- 1
    xbari <- apply(mat[cluster==i,,drop=FALSE],2, function(x) mean(x, na.rm=TRUE))
    for (j in sort(unique(clusterDel))) {
      ## measure AD
      clusi <- obsNum[cluster==i]
      clusdelj <- obsNum[clusterDel==j]
      cl <- length(clusi)*length(clusdelj)
      if (cl>0) dij[ii,jj] <- mean(matDist[clusi,clusdelj], na.rm=TRUE)
      ##      if (cl>0) dij[ii,jj] <- mean(as.matrix(Dist)[clusi,clusdelj])
      ## measure ADM
      xbarj <- apply(mat[clusterDel==j,,drop=FALSE],2, function(x) mean(x, na.rm=TRUE))
      diff <- xbari-xbarj
      if(length(diff)>0) {
        if(any.na) {
          diff <- diff[!is.na(diff)]
          dij2[ii,jj] <- sqrt(mean(diff^2))
        } else {
          dij2[ii,jj] <- sqrt(sum(diff^2))
        }
      } else {
        dij2[ii,jj] <- 0
      }
      jj <- jj+1
    }
    ii <- ii+1
  }
  rs <- matrix(rowSums(overlap),nrow=nrow(overlap),ncol=ncol(overlap),byrow=FALSE)
  cs <- matrix(colSums(overlap),nrow=nrow(overlap),ncol=ncol(overlap),byrow=TRUE)
  stabmeas["APN"] <- 1-sum(overlap^2/rs)/sum(overlap)
  stabmeas["AD"] <- sum(overlap*dij)/nrow(mat)
  stabmeas["ADM"] <- sum(overlap*dij2)/nrow(mat)
  xbar <- tapply(mat[,del],clusterDel, function(x) mean(x, na.rm=TRUE))
  stabmeas["FOM"] <- sqrt(mean((mat[,del]-xbar[as.character(clusterDel)])^2, na.rm=TRUE))/sqrt((nrow(mat)-nc1)/nrow(mat))
  return(stabmeas)
}



####################################################################
## Internal Validation Functions
####################################################################
## Silhouette
## Connectivity
## Dunn index
####################################################################



########################################################################
## Silhouette
## Value in [-1,1] to be maximized
## NOTE! cluster library also has 'silhouette' function!
## Probably use that instead
#####################################################################


mysilhouette <- function(distance=NULL, clusters, Data=NULL, method="euclidean"){

  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  dista <- apply(distance,2,function(x) tapply(x, clusters, function(x) na.rm=TRUE))
  nc <- ncol(dista); nr <- nrow(dista);
  a <- dista[matrix(c(clusters,1:nc),ncol=2,nrow=nc)]
  distb <- matrix(dista[-(clusters+(0:(nc-1))*nr)],ncol=nc,nrow=(nr-1))
  b <- apply(distb,2, function(x) min(x, na.rm=TRUE))
  s <- (b-a)/pmax(a,b)
  return(mean(s))
}



##########################################################################
## Connectivity
## Value in [0,infty] to be *minimized*
#####################################################################

connectivity <- function(distance=NULL, clusters, Data=NULL, neighbSize=10, method="euclidean"){

  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  nearest <- apply(distance,2,function(x) sort(x,ind=TRUE)$ix[2:(neighbSize+1)])
  nr <- nrow(nearest);nc <- ncol(nearest)
  same <- matrix(clusters,nrow=nr,ncol=nc,byrow=TRUE)!=matrix(clusters[nearest],nrow=nr,ncol=nc)
  conn <- sum(same*matrix(1/1:neighbSize,nrow=nr,ncol=nc))
  return(conn)
}



##########################################################################
## Dunn index
## Value in [0,infty], to be maximized
#####################################################################

dunn <- function(distance=NULL, clusters, Data=NULL, method="euclidean"){

  if (is.null(distance) & is.null(Data)) stop("One of 'distance' or 'Data' is required")
  if (is.null(distance)) distance <- as.matrix(dist(Data, method=method))
  if (class(distance)=="dist") distance <- as.matrix(distance)
  nc <- max(clusters)
  interClust <- matrix(NA, nc, nc)
  intraClust <- rep(NA, nc)

  for (i in 1:nc) {
    c1 <- which(clusters==i)
    for (j in i:nc) {
      if (j==i) intraClust[i] <- max(distance[c1,c1])
      if (j>i) {
        c2 <- which(clusters==j)
        interClust[i,j] <- min(distance[c1,c2])
      }
    }
  }
  dunn <- min(interClust,na.rm=TRUE)/max(intraClust)
  return(dunn)
}

#####################################################################
## Biological validation functions
####################################################################
## BHI (homogeneity index)
## BSI (stability index)
## Requires annotate and GO packages
#####################################################################

## read in csv file for biological annotation
readAnnotationFile <- function(filename) {
  if(!file.exists(filename))
    stop(paste(filename, "doesn't exist"))
  infile <- scan(filename, what="character")
  res <- lapply(infile, function(x) strsplit(x, ",")[[1]][-1])
  ## remove blank entries
  res <- lapply(res, function(x) x[!x%in%""])
  names(res) <- sapply(infile, function(x) strsplit(x, ",")[[1]][1])
  return(res)
}
## NOTE: returned value will be a LIST ...
## Is ok, since converted automatically to matrix in clValid function ...

## Convert annotation list to annotation TF matrix
annotationListToMatrix <- function(annotation, genenames) {
  annotation.matrix <- matrix(FALSE, ncol=length(annotation), nrow=length(genenames))
  colnames(annotation.matrix) <- names(annotation)
  rownames(annotation.matrix) <- genenames
  for ( i in 1:length(annotation) )
    {
      annot <- names(annotation)[i]
      genes <- as.character(annotation[[i]][!is.na(annotation[[i]])])
      genes.common <- intersect(genenames, genes)
      annotation.matrix[genes.common, annot ] <- TRUE
    }
  return(annotation.matrix)
}


## NOTE: could return BHI VECTOR (length=#stat clusters), take mean in clValid fn.
BHI <- function(statClust,annotation,names=NULL,category="all",dropEvidence=NULL) {

  ## Case 1
  ## Biological clusters provided by user
  if(is.matrix(annotation)) {
    ## initialize BHI vector to 0s
    bhi <- numeric(length(unique(statClust)))
    names(bhi) <- unique(statClust)

    ## for each statClust
    for ( k in unique(statClust) )
      {
        Ck.bhi <- 0
        Ck.idx <- which(statClust==k) # row indices of this statClust Ck

        if ( length(Ck.idx)<2 ) next # only one gene, skip

        ## for each gene in this statClust
        for ( i in Ck.idx )
          {
            ## ... count how many other genes j in Ck share any (1 or more)
            ## of gene i's annotations:

            B <- which(annotation[i,]==TRUE) # get indices of i's annotations
            if ( length(B)==0 ) next # gene i has no annotation

            ## gene's annotations of all other genes j in statClust Ck
            annot <- annotation[Ck.idx[ Ck.idx!= i ],B]
            ## ... add number of genes with at least one shared annotation
            if ( length(B)==1 )      Ck.bhi <- Ck.bhi + sum(annot)
            else if ( length(B) >1 ) Ck.bhi <- Ck.bhi + sum(rowSums(annot)>0)
          }

        nk <- sum(rowSums(annotation[Ck.idx,])>0) # nr. of annot. feat. in Ck
        if ( nk>1 ) bhi[k] <- Ck.bhi / (nk*(nk-1))
      }
    return(mean(bhi, na.rm=TRUE))
  }


  ## Case 2
  ## Name of annotation package provided by user
  ## Gene names assumed to correspond with rownames
  ## Requires Biobase, annotation
  ## Gene names and type of id provided by user
  ##  category <- match.arg(category,c("all","BP","CC","MF"))
  ##  expr <- parse(text="require(x, character.only=TRUE)")[[1]]
  ##  expr <- do.call("substitute", list(expr, list(x=annotation)))
  if(!require(annotation, character.only=TRUE)) {
    cat(paste("package",annotation,"not found, attempting download from Bioconductor\n",
              sep=" "))
    source("http://bioconductor.org/biocLite.R")
    res <- try(biocLite(annotation))
    if(class(res)=="try-error") {
      stop(paste("attempted download of package", annotation, "failed, exiting"))
    } else {
      library(annotation, character.only=TRUE)
    }
  }
  ##           if(!require(annotation,character.only=TRUE)) {
  ##             stop(paste("package",annotation,"not found",sep=" "))
  ##           }
  require(annotate)
  ## Get GO terms associated with each probe ID
  goTerms <- getGO(names,annotation)
  if (!is.null(dropEvidence))
    goTerms <- lapply(goTerms, dropECode, dropEvidence)
  ## Find biological homogeneity for each stat cluster, then take average
  bhi <- tapply(goTerms,statClust,function(x)  matchGO(x,category))
  return(mean(bhi, na.rm=TRUE))
  ## prev: return(mean(bhi[bhi!=-9]))
} ## End BHI function


matchGO <- function(gg,category) {
  ## x[1] is x$GOID, x[3] is x$Ontology
  ## extracts the GO_ID and Ont for each GO term associated with each probe ID
  goIDs <- lapply(gg, function(a) sapply(a, function(x) x[1]))
  ont   <- lapply(gg, function(a) sapply(a, function(x) x[3]))
  switch(category,
         ## 1. Find which probes have corresponding annotation
         ## 2. Extract only subset w/appropriate annotation
         BP = {
           goBP <- sapply(ont, function(x) any(x%in%"BP"))
           goIDs <- goIDs[goBP]
         },
         CC = {
           goCC <- sapply(ont, function(x) any(x%in%"CC"))
           goIDs <- goIDs[goCC]
         },
         MF = {
           goMF <- sapply(ont, function(x) any(x%in%"MF"))
           goIDs <- goIDs[goMF]
         },
         all = {
           ## goAll <- sapply(goIDs, function(a) all(sapply(a,function(x) is.null(x))))
           ## changed 11/17/09
           goAll <- sapply(ont, function(x) any(x%in%c("BP","CC","MF")))
           goIDs <- goIDs[goAll]
         })
  ## Number of annotated probes
  ## If only one annotated probe, then skip this cluster
  n <- length(goIDs)
  if (n<2) return(NA) ## previously: return(-9)
  sum <- 0
  ## Now count overlap of annotation categories for different probes/genes
  for (i in 1:(length(goIDs)-1)) {
    for (j in (i+1):length(goIDs)) {
      switch(category,
             all =  sum <- sum + any(goIDs[[i]]%in%goIDs[[j]]),
             BP  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="BP"]%in%goIDs[[j]][ont[[j]]=="BP"]),
             CC  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="CC"]%in%goIDs[[j]][ont[[j]]=="CC"]),
             MF  =  sum <- sum + any(goIDs[[i]][ont[[i]]=="MF"]%in%goIDs[[j]][ont[[j]]=="MF"]))
    }
  }
  return(sum/(n*(n-1)))
} ## End matchGO function


BSI <- function(statClust,statClustDel,annotation,names=NULL,category="all",
                goTermFreq=0.05, dropEvidence=NULL) {

  ## Case 1
  ## Biological clusters provided by user
  if(is.matrix(annotation)) {
    nFC <- ncol(annotation)
    nAnnot <- apply(annotation, 2, sum)
    if(is.null(names(statClust))) {
      names(statClust) <- names
      names(statClustDel) <- names
    }
    bsi <- numeric(nFC)
    overlap <- xtabs(~statClust + statClustDel)
    rsums <- rowSums(overlap)
    for (FCk in 1:nFC)
      {
        osum <- 0
        g.idx <- which(annotation[,FCk]) ## find which genes are annotated
        if (length(g.idx)<2) next ## only one gene, skip
        for (gx in g.idx) {
          for (gy in g.idx) {
            if (gx != gy) {
              i <- statClust[gx]
              j <- statClustDel[gy]
              osum <- osum + overlap[i,j]/rsums[i]
            }
          }
        }
        bsi[FCk] <- osum/(nAnnot[FCk]*(max(nAnnot[FCk]-1,1)))
      }
    return(mean(bsi, na.rm=TRUE))
  }

  ## Case 2
  ## Gene names and type of id provided by user
  ## For option 2 of BSI need to determine how many GO terms to use
  ## Cutoff determined by goTermFreq (cuts all GO terms w/ < 5% freq in Data set)
  ##  category <- match.arg(category,c("all","BP","CC","MF"))

  tab <- xtabs(~statClust + statClustDel)
  rs <- rowSums(tab)
  n <- length(statClust)


  if(!require(annotation,character.only=TRUE)) {
    cat(paste("package",annotation,"not found, attempting download from Bioconductor\n",
              sep=" "))
    source("http://bioconductor.org/biocLite.R")
    try(biocLite(annotation))
  }
  ##           if(!require(annotation,character.only=TRUE)) {
  ##             stop(paste("package",annotation,"not found",sep=" "))
  ##           }
  goTerms <- getGO(names,annotation)
  if (!is.null(dropEvidence))
    goTerms <- lapply(goTerms, dropECode, dropEvidence)

  ## Things to do
  ## 1. extract all relevant GO terms for each gene (bp,cc, etc)
  ## 2. get freq table of terms - keep most frequent
  ## 3. create indicator matrix for each term indicating which genes have that term

  switch(category,
         ## x[1] is x$GOID, x[3] is x$Ontology
         BP = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="BP"])),
         CC = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="CC"])),
         MF = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1][x[3]=="MF"])),
         all = goIDs <- sapply(goTerms, function(a) sapply(a, function(x) x[1])))

  goTab <- table(unlist(goIDs))
  keepTerms <- names(goTab)[goTab>floor(n*goTermFreq)]
  termMat <- matrix(0,ncol=length(keepTerms),nrow=n)
  ## 09/27/09 - fixed issue by using unlist(x) in sapply
  for (i in 1:length(keepTerms)) {
    termMat[,i] <- sapply(goIDs, function(x) keepTerms[i]%in%unlist(x))
  }

  bsi <- apply(termMat,2, function(a) {
    out <- outer(statClust[as.logical(a)],statClustDel[as.logical(a)], function(x,y) tab[cbind(x,y)]/rs[x])
    (sum(out) - sum(diag(out)))/(sum(a)*(ifelse(sum(a)>1,sum(a)-1,1)))
  })

  return(mean(bsi, na.rm=TRUE))

}


#####################################################################
## Functions for clustering - SOTA
#####################################################################


sota.init <- function(data){
  nodes <- matrix(0, nrow(data)*2, 3+ncol(data))
  if(is.null(colnames(data)))
    colnames(data) <- paste("V", 1:ncol(data))
  colnames(nodes) <- c("ID", "anc", "cell", colnames(data))
  nodes[,"ID"]=1:(nrow(data)*2)

  nodes[1,] <- c(1, 0, 0, apply(data,2, function(x) mean(x, na.rm=TRUE)))
  nodes[2,] <- c(2, 1, 1, nodes[1,][-c(1,2,3)])
  nodes[3,] <- c(3, 1, 1, nodes[1,][-c(1,2,3)])
  return(nodes)
}

dist.fn <- function(input, profile, distance){
  if(distance=="correlation")
    return(1-cor(input,profile, use="pairwise.complete.obs"))
  else
    return(sqrt(sum((input-profile)^2)))
}

cl.ID <- function(clust, old.cl, new.cl){
  for(i in 1:length(clust))
    clust[i] <- new.cl[which(old.cl==clust[i])]
  clust
}

getResource <- function(data, tree, clust, distance, pr){
  dist <- rep(0, length(clust))
  resource <- rep(0, max(clust))

  for(i in unique(clust)){
    temp <- data[clust==i,]
    if(is.vector(temp))
      temp <- matrix(temp, nrow=1, ncol=ncol(data))
    if(distance=="correlation")
      resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr],
                                distance=distance))
    else
      resource[i] <- mean(apply(temp, 1, dist.fn, profile=tree[i,pr],
                                distance=distance))}
  resource
}

getCells <- function(tree, neighb.level, n){
  or.n <- n
  cells <- c(n-1,n)
  for(i in 1:(neighb.level+1)){
    n  <- tree[n, "anc"]
    if(n==1)
      break
  }
  for(j in 2:(or.n-2)){
    z <- j
    if(tree[j,"cell"]!=1)
      next
    while(z > 0){
      z <- tree[z, "anc"]
      if(z==n){
        cells <- c(cells, j)
        break}
    }
  }
  return(tree[cells,])
}


sota <- function(data, maxCycles, maxEpochs=1000, distance="euclidean",
                 wcell=.01, pcell=.005, scell=.001, delta=.0001, neighb.level=0,
                 maxDiversity = .9, unrest.growth=TRUE, ...){
  tree <- sota.init(data)
  pr <- 4:ncol(tree)
  n <- 3
  genes<- 1:nrow(data)
  clust <- rep(1, length(genes))
  Node.Split <- 1

  Res.V <- getResource(data, tree, clust, distance, pr)
  if(distance=="correlation")
    Res.V <- 1-Res.V
  diversity <- Res.V

  for(k in 1:maxCycles){                                #loop for the Cycles
    trainNode <- Node.Split
    trainSamp <- genes[clust==trainNode]
    curr.err <- 1e10
    ep <- 1

    while(ep <= maxEpochs){	      			#loop for the Epochs
      last.err <- 0
      left.ctr <- right.ctr <- 0
      left.d <- right.d <- 0
      for(i in trainSamp){
        cells <- tree[c(n-1,n),]
        dist <- rep(0, nrow(cells))
        for(j in 1:2)
          dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)

        or <- which.min(dist)
        if(or==1)
          left.ctr <- left.ctr + 1
        else
          right.ctr<- right.ctr + 1

        closest <- cells[or,1]
        sis <- ifelse(closest%%2==0,closest+1,closest-1)
        sis.is.cell <- ifelse(tree[sis,"cell"]== 1, 1, 0)

        ##   updating the cell and its neighbourhood
        if(sis.is.cell==1){
          parent <- tree[closest, "anc"]
          tree[closest, pr] <- tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
          tree[sis, pr] <- tree[sis, pr]+scell*(data[i,]-tree[sis,pr])
          tree[parent, pr] <- tree[parent, pr]+pcell*(data[i,]-tree[parent, pr])
        }
        else
          {
            tree[closest, pr] <- tree[closest, pr]+wcell*(data[i,]-tree[closest, pr])
          }
      }
      cells <- tree[c(n-1,n),]
      for(i in trainSamp){
        for(j in 1:2)
          dist[j] <- dist.fn(data[i,], cells[j,pr], distance=distance)
        last.err <- last.err+min(dist)}
      last.err <- last.err/length(trainSamp)

      if(ifelse(last.err==0, 0, abs((curr.err-last.err)/last.err)) < delta
         && left.ctr !=0 && right.ctr !=0)
        break
      ep <- ep + 1
      curr.err <- last.err
    }
    clust <- assignGenes(data, trainSamp, clust, tree, n, distance, pr, neighb.level)
    Res.V <- getResource(data, tree, clust, distance, pr)
    if(distance=="correlation")
      Res.V <- 1-Res.V

    tempRes <- Res.V
    tempRes[tempRes == 0] <- diversity[tempRes==0]
    diversity <- tempRes

    if(k==maxCycles || (max(Res.V) < maxDiversity & unrest.growth==FALSE))
      break  ## do not split the cell
    newCells <- splitNode(Res.V, tree, n)
    tree <- newCells$tree
    n <- newCells$n
    Node.Split <- newCells$toSplit
  }

  tree <- trainLeaves(data, tree, clust, pr, wcell, distance, n, delta)
  Res.V <- getResource(data, tree, clust, distance, pr)
  Res.V <- Res.V[Res.V!=0]
  if(distance=="correlation")
    Res.V <- 1-Res.V

  diversity[(length(diversity)-length(Res.V)+1):length(diversity)] <- Res.V

  treel <- tree[tree[,"cell"]==1,]
  old.cl <- treel[,1]
  treel[,1] <- 1:nrow(treel)
  old.clust <- clust

  clust <- cl.ID(old.clust, old.cl, 1:nrow(treel))
  totals <- table(clust)

  out <- list(data=data, c.tree=cbind(tree[1:n,],Diversity=diversity), tree=treel, clust=clust,
              totals=totals, dist=distance, diversity=Res.V)

  class(out) <- "sota"
  return(out)
}



trainLeaves <- function(data, tree, clust, pr, wcell, distance, n, delta){
  nc <- ncol(data)
  for(i in 1:n){
    if(!is.element(i, clust))
      next
    temp <- matrix(data[clust==i,], ncol=nc)
    converged <- FALSE
    init.err <- getCellResource(temp, tree[i,pr], distance)
    while(!converged){
      for(j in 1:nrow(temp))
        tree[i, pr] <- tree[i, pr]+wcell*(temp[j,]-tree[i, pr])

      last.err <- getCellResource(temp, tree[i,pr], distance)
      converged <- ifelse(abs((last.err-init.err)/last.err) < delta, TRUE, FALSE)
      init.err <- last.err
    }
  }
  return(tree)
}


assignGenes <- function(data, Sample, clust, tree, n, distance, pr, neighb.level){
  if(neighb.level==0)
    cells <- tree[c(n-1,n),]
  else
    cells <- getCells(tree, neighb.level, n)

  for(i in Sample){
    dist <- rep(0, nrow(cells))
    for(j in 1:nrow(cells))
      dist[j] <- dist.fn(data[i,], cells[j,pr], distance)
    or <- which.min(dist)
    closest <- cells[or,1]
    clust[i] <- closest
  }
  clust
}

splitNode <- function(Res.V, tree, n){
  maxheter <- which.max(Res.V)
  cl.to.split <- tree[maxheter,1]
  tree[n<-n+1,-1] <- tree[cl.to.split,-1]
  tree[n, "anc"] <- cl.to.split
  tree[n<-n+1,-1] <- tree[cl.to.split,-1]
  tree[n, "anc"] <- cl.to.split
  tree[cl.to.split, "cell"] <- 0
  return(list(tree=tree, n=n, toSplit=cl.to.split))
}



getCellResource <- function(temp, profile, distance){
  if(distance=="correlation")
    resource <- mean(apply(temp, 1, dist.fn, profile,
                           distance=distance))
  else(distance=="euclidean")
  resource <- mean(apply(temp, 1, dist.fn, profile,
                         distance=distance))
  resource
}

print.sota <- function(x, ...){
  results <- as.matrix(cbind(as.numeric(names(x$totals)), as.numeric(x$totals),
                             x$diversity))  ## changed
  colnames(results) <- c("ID","Size", "Diversity")
  rownames(results) <- rep("", nrow(results))
  cat("\nClusters:\n")
  print(results, ...)
  cat("\nCentroids:\n")

  print(format(data.frame(x$tree[,-c(1:3)])), ...)
  cat("\n")
  cat(c("Distance: ", x$dist, "\n"))
  invisible(x)
}


plot.sota <- function(x, cl=0, ...){

  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  if(cl!=0)
    par(mfrow=c(1,1)) else
  {
    pdim <- c(0,0)
    for(i in 1:100){
      j <- i
      if(length(x$totals) > i*j)
        j <- j+1
      else{
        pdim <- c(i,j)
        break}
      if(length(x$totals) > i*j)
        i <- i+1
      else{
        pdim <- c(i,j)
        break}
    }
    par(mfrow=pdim)
  }

  ylim = c(min(x$data), max(x$data))
  pr <- 4:ncol(x$tree)
  if(cl==0)
    cl.to.print <- 1:length(table(x$clust)) else  ## changed
  cl.to.print <- cl
  cl.id <- sort(unique(x$clust))  ## changed

  for(i in cl.to.print){
    plot(1:ncol(x$data), x$tree[i, pr], col="red", type="l",
         ylim=ylim, xlab=paste("Cluster ",i), ylab="Expr. Level", ...)
    legend("topleft", legend=paste(x$totals[i], " Genes"), cex=.7,
           text.col="navy", bty="n")
    cl <- x$data[x$clust==cl.id[i],]  ## changed
    if(is.vector(cl))
      cl <- matrix(cl, nrow=1)
    for(j in 1:x$totals[i])
      lines(1:ncol(x$data), cl[j,], col="grey")
    lines(1:ncol(x$data), x$tree[i, pr], col="red", ...)

  }
}


##################################################################################
## Functions for Rank Aggregation
## 03/08/2009
## NOTE: Will encorporate this into RankAggreg instead
##################################################################################

getRanksWeights <- function(clVObj, measures=measNames(clVObj), nClust=nClusters(clVObj),
                            clAlgs=clusterMethods(clVObj)) {
  measures <- match.arg(measures,measNames(clVObj),several.ok=TRUE)
  nClust <- as.character(nClust)
  nClust <- match.arg(nClust,nClusters(clVObj),several.ok=TRUE)
  nClust <- as.character(nClust)
  clAlgs <- match.arg(clAlgs,clusterMethods(clVObj),several.ok=TRUE)
  meas <- clVObj@measures[measures, nClust, clAlgs, drop=FALSE]
  ranks <- matrix(0, dim(meas)[1], dim(meas)[2]*dim(meas)[3])
  colnames(ranks) <- 1:ncol(ranks)
  rownames(ranks) <- rownames(meas)
  weights <- ranks

  tmp <- expand.grid(dimnames(meas)[[2]], dimnames(meas)[[3]])
  algs <- paste(tmp[,2], tmp[,1], sep="-")
  ## measures <- rownames(meas)

  for(i in 1:nrow(meas)){
    if(measures[i] %in% c("Dunn", "Silhouette"))
      incr <- TRUE
    else
      incr <- FALSE
    sorted <- sort(meas[i,,], ind=TRUE, decreasing=incr)
    ranks[i,] <- algs[sorted$ix]
    weights[i,] <- sorted$x
  }
  list(ranks=ranks, weights=weights)
}
