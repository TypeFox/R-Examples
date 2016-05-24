

## crimeClust_hier
##==============================================================================
#' Agglomerative Hierarchical Crime Series Clustering
#'
#' Run hierarchical clustering on a set of crimes using the log Bayes Factor as 
#' the similarity metric.
##  Inputs:
#'  @param crimedata data.frame of crime incidents. Must contain a column named
#'   \code{crimeID}.
#'  @param varlist a list of the variable names (columns of \code{crimedata}) 
#'   used to create evidence variables with \code{\link{compareCrimes}}. 
#'  @param estimateBF function to estimate the log bayes factor from evidence 
#'   variables
#'  @param linkage the type of linkage for hierarchical clustering
#'    \itemize{
#'      \item \dQuote{average} uses the average bayes factor
#'      \item \dQuote{single} uses the largest bayes factor (most similar)
#'      \item \dQuote{complete} uses the smallest bayes factor (least similar)
#'    }
#'  @param \ldots other arguments passed to \code{\link{compareCrimes}}
#'  @details This function first compares all crime pairs using \code{\link{compareCrimes}},
#'    then uses \code{estimateBF} to estimate the log Bayes factor for every pair.
#'    Next, it passes this information into \code{\link{hclust}} to carry out the 
#'    agglomerative hierarchical clustering. Because \code{\link{hclust}} requires 
#'    a dissimilarity, this uses the negative log Bayes factor. 
#'    
#'    The input \code{varlist} is a list with elements named: crimeID, spatial, 
#'    temporal, categorical, and numerical. Each element should be a vector of 
#'    the column names of \code{crimedata} corresponding to that feature. See 
#'    \code{\link{compareCrimes}} for more details. 
##  Outputs:
#'  @return An object of class \code{hclust} (from \code{\link{hclust}}). 
#'  @seealso \code{\link{clusterPath}, \link{plot_hcc}}
#'  @examples
#'  data(crimes)
#'  #- cluster the first 10 crime incidents
#'  crimedata = crimes[1:10,]
#'  varlist = list(spatial = c("X", "Y"), temporal = c("DT.FROM","DT.TO"), 
#'      categorical = c("MO1",  "MO2", "MO3"))
#'  estimateBF <- function(X) rnorm(NROW(X))   # random estimation of log Bayes Factor
#'  HC = crimeClust_hier(crimedata,varlist,estimateBF)
#'  plot_hcc(HC,yticks=-2:2)
#'  
#'  # See vignette: "Crime Series Identification and Clustering" for more examples.
#'  @references
#'  Porter, M. D. (2014). A Statistical Approach to Crime Linkage.
#'    \emph{arXiv preprint arXiv:1410.2285.}.
#'  \url{http://arxiv.org/abs/1410.2285}
#'  @export
##  #ToDO: add crime.labels = ifknown, crimegroups  
##==============================================================================
crimeClust_hier <- function(crimedata,varlist,estimateBF, 
                        linkage = c('average','single','complete'),...){
  linkage = match.arg(linkage)
  crimeIDs = unique(as.character(crimedata$crimeID))  
  allPairs = t(combn(crimeIDs,2))
  A = compareCrimes(allPairs,crimedata,varlist,...)
  bf = estimateBF(A)  
  #-- Perform aglomerative hierarchcical clustering
  d2 = as.numeric(-bf); class(d2) = 'dist'
  attr(d2,"method") = "log Bayes Factor"
  attr(d2,"Labels") = crimeIDs 
  attr(d2,"Size") = length(crimeIDs)
  offset = ceiling(max(bf))
  d2 = d2 + offset
  hc = switch(linkage,
                      'single'   = hclust(d2,method='single'),
                      'complete' = hclust(d2,method='complete'),
                      'average'  = hclust(d2,method='average'))
  hc$offset = offset
return(hc)
}




## plot_hcc
##==============================================================================
#' Plot a hierarchical crime clustering object
#'
#' Similar to \code{\link{plot.dendrogram}}. 
##  Inputs:
#'  @param tree an object produced from \code{\link{crimeClust_hier}}
#'  @param yticks the location of the tick marks for log Bayes factors
#'  @param hang the hang argument of \code{\link{as.dendrogram}}
#'  @param \ldots other arguments passed to \code{\link{plot.dendrogram}}
#'  @details This function creates a dendrogram object and then plots it. It
#'    corrects the y-axis to give the proper values and adds the number of clusters
#'    if the tree were cut at a particular log Bayes factor.
##  Outputs:
#'  @return  A dendrogram
#'  @examples
#'  # See vignette: "Crime Series Identification and Clustering" for usage.
#'  @seealso \code{\link{crimeClust_hier}}
#'  @export
##==============================================================================
plot_hcc <- function(tree,yticks=seq(-2,8,by=2),hang=-1,...){
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mar=c(3.5,4,3.15,4)+.1)
  hcd = as.dendrogram(tree,hang=hang)  
  plot(hcd,yaxt="n",...) 
  offset = tree$offset
  labs = -yticks
  abline(h=labs+offset,col="grey80",lty=3)
  axis(2,at=labs+offset,labels=-labs,cex=.6,las=1)
  nClusters = sapply(labs+offset,function (h) length(unique(cutree(tree,h=h))))
  axis(4,at=labs+offset,labels=nClusters,las=1)
  title(ylab='log Bayes factor')
  mtext('number of clusters',side=4,line=3)
}


## clusterPath
##==============================================================================
#' Follows path of one crime up a dendrogram
#'
#'  The sequence of groups that a crime belongs to. 
##  Inputs:
#'  @param crimeID the crime ID for a crime used in hierarchical clustering
#'  @param tree an object produced from \code{\link{crimeClust_hier}}
#'  @details Agglomerative hierarchical clustering form clusters by sequentially
#'    merging the most similar groups at each iteration. This function is designed
#'    to help trace the sequence of groups an individual crime is a member of. And
#'    it shows at what score (log Bayes factor) the merging occurred.
##  Outputs:
#'  @return  data.frame of the additional crimes and the log Bayes factor at each
#'    merge.
#'  @examples
#'  # See vignette: "Crime Series Identification and Clustering" for usage.    
#'  @seealso \code{\link{crimeClust_hier}, \link{plot_hcc}}
#'  @export
##==============================================================================
clusterPath <- function(crimeID,tree){ 
  crimeID = as.character(crimeID)
  bf = -tree$height + tree$offset          # actual log bayes factor 
  
  #-- Find rows of tree$merge corresponding to event of interest
  findVal <- function(id,X=tree$merge) which(X[,1]==id | X[,2]==id)
  id = which(tree$labels == crimeID)       # label index
  ind = i = findVal(-id)
  repeat{
    #i = which(tree$merge[,1]==i | tree$merge[,2]==i)
    i = findVal(i)
    ind = c(ind,i)
    if(length(i)==0) break
  }    

  #-- List of element at each level of tree
  getElements <- function(x){
    if(x[2] < 0)    tree$labels[-x]
    else if(x[1]<0) c(treeList[[x[2]]],tree$labels[-x[1]])
    else            unlist(treeList[x])
  }
  treeList = list()
  for(i in 1:length(tree$height)){
    treeList[[i]] = getElements(tree$merge[i,])
  }
  tL = treeList[ind]                  # only elements that include crimeID

  #-- Find which new elements are added at each merge
  new.crimes = setdiff(tL[[1]],crimeID)
  for(i in 2:length(tL)){
    d = setdiff(tL[[i]],tL[[i-1]])
    new.crimes[i] = paste(sort(d),collapse=', ') 
  }
  
  DF = data.frame(logBF=bf[ind],crimes=new.crimes) 
  attr(DF,"crimeID") = crimeID
return(DF)
}
