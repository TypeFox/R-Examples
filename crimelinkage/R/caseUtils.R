
################################################################################
####  caseUtils.R
####  Various functions to generate data for case Linkage
##
## offenderTable - data.frame of solved crimes matching crimeID to offenderID
##     has columns named: offenderID, crimeID
## crimeID - the ID for a crime (class: character)
## offenderID - the ID for an offender (class: character)
################################################################################

#-- Required R Packages
#require(igraph)  # for makeUnlinked()



## getCriminals
##==============================================================================
#' Lookup the offenders responsible for a set of solved crimes
#'
#' Generates the IDs of criminals responsible for a set of solved crimes using
#' the information in \code{offenderTable}.
##  Inputs:
#'  @param crimeID crimeID(s) of solved crimes.
#'  @param offenderTable offender table that indicates the offender(s) responsible 
#'    for solved crimes. \code{offenderTable} must have columns named: 
#'    \code{offenderID} and \code{crimeID}.
##  Outputs:
#'  @return Vector of offenderIDs responsible for crimes labeled \code{crimeID}.
#'  @seealso \code{\link{getCrimeSeries}}
#'  @examples
#'  
#'  data(offenders)
#'  
#'  getCriminals("C:1",offenders)
#'  
#'  getCriminals("C:78",offenders)                      # shows co-offenders
#'  
#'  getCriminals(c("C:26","C:78","85","110"),offenders) # all offenders from a crime series
#'  @export
##==============================================================================
getCriminals <- function(crimeID,offenderTable){
  crimeID = as.character(crimeID)
  ind = which(offenderTable$crimeID %in% crimeID)
  offenderID = unique(as.character(offenderTable$offenderID[ind]))  # offenders/suspects involved in crimes
  return(sort(offenderID))
}

## getCrimeSeries
##==============================================================================
#' Generate a list of offenders and their associated crime series.
##  Inputs:
#'  @param offenderID vector of offender IDs
#'  @param offenderTable offender table that indicates the offender(s) responsible 
#'    for solved crimes. \code{offenderTable} must have columns named: 
#'    \code{offenderID} and \code{crimeID}.
#'  @param restrict if vector of \code{crimeID}, then only include those crimeIDs
#'    in \code{offenderTable}. If \code{NULL}, then return all crimes for offender.
#'  @param show.pb (logical) should a progress bar be displayed
##  Outputs:
#'  @return List of offenders with their associated crime series.
#'  @seealso \code{\link{makeSeriesData}}, \code{\link{getCriminals}}, 
#'    \code{\link{getCrimes}}
#'  @examples 
#'  
#'  data(offenders)
#'  
#'  getCrimeSeries("O:40",offenders)
#'  getCrimeSeries(c("O:40","O:3"),offenders)  # list of crime series from multiple offenders
#'  @export
##==============================================================================
getCrimeSeries <- function(offenderID,offenderTable,restrict=NULL,show.pb=FALSE){
  offenderID = unique(as.character(offenderID))
  n = length(offenderID)
  CS = vector('list',n)
  valid.crimeID = TRUE
  if(!is.null(restrict)){           # Pre-compute valid crimeIDs
    valid.crimeID = (offenderTable$crimeID %in% as.character(restrict))
  }
  valid.offenderID = (offenderTable$offenderID %in% offenderID)   # Pre-compute valid offenderID number
  valid = valid.crimeID & valid.offenderID
  offenderTable = unique(offenderTable[valid,])  # only use unique valid offenderTable to speed calculations
  if(show.pb) pb = txtProgressBar(style=3,max=n)
  for(i in 1:n){
    oid = offenderID[i]
    ind = which(offenderTable$offenderID %in% oid)
    cid = as.character(offenderTable$crimeID[ind])
    CS[[i]] = list(offenderID = oid, crimeID = cid)
    if(show.pb) setTxtProgressBar(pb,i)
  }
  if(show.pb) close(pb)
  if(i==1) CS = CS[[1]]
  return(CS)
}


## getCrimes
##==============================================================================
#' Generate a list of crimes for a specific offender
##  Inputs:
#'  @param offenderID an offender ID that is in \code{offenderTable}
#'  @param offenderTable offender table that indicates the offender(s) responsible 
#'    for solved crimes. \code{offenderTable} must have columns named: 
#'    \code{offenderID} and \code{crimeID}.
#'  @param crimedata data.frame of crime incident data. \code{crimedata} must be
#'    a data.frame with a column named: \code{crimeID} 
##  Outputs:
#'  @return The subset of crimes in \code{crimedata} that are attributable to 
#'    the offender named \code{offenderID}
#'  @seealso \code{\link{getCrimeSeries}}
#'  @examples
#'  data(crimes)
#'  data(offenders)
#'  
#'  getCrimes("O:40",crimes,offenders)
#'  @export
##==============================================================================
getCrimes <- function(offenderID,crimedata,offenderTable){
  cid = offenderTable$crimeID[offenderTable$offenderID %in% offenderID]
  crimes = crimedata[crimedata$crimeID %in% cid,]
  return(crimes)
}



## makeSeriesData
##==============================================================================
#'  Make crime series data
#'
#'  Creates a data frame with index to crimedata and offender information. It is 
#'  used to generate the linkage data. 
##  Inputs:
#'  @param crimedata data.frame of crime incident data. \code{crimedata} must have
#'    columns named: \code{crimeID}, \code{DT.FROM}, and \code{DT.TO}. Note: if
#'    crime timing is known exactly (uncensored) than only \code{DT.FROM} is 
#'    required.
#'  @param offenderTable offender table that indicates the offender(s) responsible 
#'    for solved crimes. \code{offenderTable} must have columns named: 
#'    \code{offenderID} and \code{crimeID}.
#'  @param time the event time to be returned: 'midpoint', 'earliest', or
#'    'latest'  
##  Outputs:
#'  @return data frame representation of the crime series present in the 
#'    \code{crimedata}. It includes the crime ID (\code{crimeID}), index of that
#'    crimeID in the original \code{crimedata} (\code{Index}), the crime series 
#'    ID (\code{CS}) corresponding to each \code{offenderID}, and the event time 
#'    (\code{TIME}).
#'  @details The creates a crimeseries data object that is required for creating
#'    linkage data. It creates a crime series ID (\code{CS}) for every 
#'    offender. Because of co-offending, a single crime (\code{crimeID}) can 
#'    belong to multiple crime series. 
#'  @seealso \code{\link{getCrimeSeries}}
#'  @examples  
#'  data(crimes)
#'  data(offenders)
#'  
#'  seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
#'  head(seriesData)
#'  
#'  nCrimes = table(seriesData$offenderID)  # length of each crime series
#'  table(nCrimes)                  # distribution of crime series length
#'  mean(nCrimes>1)                 # proportion of offenders with multiple crimes
#'  
#'  nCO = table(seriesData$crimeID) # number of co-offenders per crime
#'  table(nCO)                      # distribution of number of co-offenders
#'  mean(nCO>1)                     # proportion of crimes with multiple co-offenders
#'  @export
#
##  Note: could use dplyr
##  require(dplyr)
##  seriesData = distinct(offenderTable) %>% 
##               inner_join(crimedata,by="crimeID") %>%
##               mutate(TIME=DT.FROM + (DT.TO-DT.FROM)/2) %>% 
##               select(crimeID,offenderID,TIME)
##==============================================================================
makeSeriesData <- function(crimedata,offenderTable,time=c("midpoint","earliest","latest")){
  cid = unique(crimedata$crimeID)
  oid = getCriminals(cid,offenderTable)   # offenders responsible for the crimes
  CS = getCrimeSeries(oid,offenderTable,restrict=cid,show.pb=FALSE)
  nCS = length(CS)
  nCrimes = sapply(CS,function(x) length(x$crimeID))
  # tabCrimes = table(nCrimes)   # number of crimes per series
  a = unlist(sapply(CS,'[','crimeID'))  # All crimes by offenders in oid crimeseries'
  b = match(a,crimedata$crimeID)      # index to Crimes
  time = match.arg(time)
  if(is.null(crimedata$DT.TO)) crimedata$DT.TO = crimedata$DT.FROM  # for uncensored data
  TIME = switch(time,
                 midpoint = with(crimedata[b,],DT.FROM + (DT.TO-DT.FROM)/2),
                 earliest = with(crimedata[b,],DT.FROM),
                 latest   = with(crimedata[b,],DT.TO)  
                )
  seriesData = data.frame(crimeID=a, Index=b, 
                          CS=rep(1:nCS,times=nCrimes),
                          offenderID=rep(oid,times=nCrimes),
                          TIME,
                          stringsAsFactors=FALSE)
return(seriesData)
}







##  makeGroups
##==============================================================================
#'  Generates crime groups from crime series data
#'
#'  This function generates crime groups that are useful for making unlinked pairs 
#'  and for agglomerative linkage.
##  Inputs:
#'  @param X crime series data (generated from \code{\link{makeSeriesData}}) 
#'    with offender ID (\code{offenderID}),
#'    crime ID (\code{crimeID}), and the event datetime (\code{TIME})
#'  @param method Method=1 (default) forms groups by finding
#'    the maximal connected offender subgraph. Method=2 forms groups from the 
#'    unique group of co-offenders. Method=3 forms from groups from offenderIDs
##  Outputs:
#'  @return vector of crime group labels 
##  Notes:
#'  @details Method=1 forms groups by finding the maximal connected offender 
#'    subgraph. So if two offenders have ever co-offended, then all of their crimes
#'    are assigned to the same group. Method=2 forms groups from the unique group 
#'    of co-offenders. So for two offenders who co-offended, all the co-offending 
#'    crimes are in one group and any crimes committed individually or with other 
#'    offenders are assigned to another group. Method=3 forms groups from the 
#'    offender(s) responsible. So a crime that is committed by multiple people 
#'    will be assigned to multiple groups.
#'  @examples 
#'   
#'  data(crimes)
#'  data(offenders)
#'  seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
#'  groups = makeGroups(seriesData,method=1)
#'  head(groups,10)
#'  
#'  @export
##  igraph package is needed for getting crime groups for Method==1
##==============================================================================
makeGroups <- function(X,method=1){ 
  
  if(method==1){ # Maximal connected offender subgraph
    # Group formed from maximal connected offender subgraphs
    #  Problem: Some groups will be large, and possibly have dissimilar crimes
    #  Benefit: There is no overlap in criminals or crimes between groups

    # require(igraph)
    #-- Make edgelist of connected offenders
    pairwise <- function(A) if(length(A)>1) t(combn(A,2))
    Y = tapply(X$offenderID,X$crimeID,unique)        # list of offender IDs for each crime
    EL = do.call(rbind,sapply(Y,pairwise))  # edgelist for connected offenders
    #-- Get connected sets (graph clusters)
    G = unique(X$offenderID)        # unique offenders labels
    Graph = igraph::simplify(igraph::graph.data.frame(EL,directed=FALSE,
                                                      vertices=data.frame(G)))
    Graph.clusters = igraph::clusters(Graph)
    C = Graph.clusters$membership
    nCG = length(unique(C))  # number of crime groups
    #-- make crime group data
    CG = C[match(X$offenderID,igraph::V(Graph)$name)]   # crime group for each row in X
  }
  
  #------------------------------
  if(method==2){ # Unique group of co-offenders
    # Group formed from crimes committed by same combination of co-offenders
    #  Problem: some groups will contain same offenders and thus should be similar
    #  Benefit: groups of different offenders might have different behavior

    ID = tapply(X$offenderID,X$crimeID,function(x) paste(sort(unique(x)),collapse=", "))  
    g = unique(data.frame(crimeID=names(ID),group=ID,stringsAsFactors=FALSE))
    g$CG = match(g$group,unique(ID))
    # g = seriesData %>% group_by(crimeID) %>% 
    #  summarize(group=paste(sort(unique(offenderID)),collapse=', '))
    # g$CG = match(b$group,unique(b$group))
    CG = g$CG[match(X$crimeID,g$crimeID)]
  }
  
  #------------------------------
  if(method==3){ # Each offender is a group
    # Group formed from each offender
    #  Problem: Crimes will belong to multiple groups
    #  Benefit: Each group is specific to a single offender
    
    CG = match(X$offenderID,unique(X$offenderID))
  } 
  CG = as.integer(CG)
return(CG)
}













##  makePairs
##==============================================================================
#'  Generates indices of linked and unlinked crime pairs (with weights)
#'
#'  These functions generate a set of crimeIDs for linked and unlinked crime pairs.
#'    Linked pairs are assigned a weight according to how many crimes are in the 
#'    crime series. For unlinked pairs, \code{m} crimes are selected from each 
#'    \emph{crime group} and pairs them with crimes in other \emph{crime groups}.  
##  Inputs:
#'  @param X crime series data (generated from \code{\link{makeSeriesData}}) 
#'    with offender ID (\code{offenderID}),
#'    crime ID (\code{crimeID}), and the event datetime (\code{TIME})
#'  @param thres the threshold (in days) of allowable time distance
#'  @param m the number of samples from each crime group (for unlinked pairs)
#'  @param show.pb (logical) should a progress bar be displayed
#'  @param seed seed for random number generation
##  Outputs:
#'  @return matrix of indices of crime pairs with weights. For \code{makePairs},
#'   The last column \code{type} indicates if the crime pair is linked or unlinked.
##  Notes:
#'  @details 
#'   \code{makePairs} is a Convenience function that calls \code{makeLinked} and
#'   \code{makeUnlinked} and combines the results. It is unlikely that the latter
#'   two functions will need to be called directly.
#'   
##  makeLinked:  
#'   For \emph{linked} crime pairs, the weights are such that each crime series
#'   contributes a total weight of no greater than 1. Specifically, the weights 
#'   are  \eqn{W_{ij} = \min \{1/N_m: V_i,V_j \in C_m \}}, 
#'   where \eqn{C_m} is the crime series for offender \eqn{m} and \eqn{N_m} is 
#'   the number of crime pairs in their series (assuming \eqn{V_i} and \eqn{V_j} 
#'   are together in at least one crime series).
##   such that each crime
##   series contributes a total weight of 1. 
#'   Due to co-offending, the sum of weights 
#'   will be smaller than the number of series with at least two crimes.
#'   
##  makeUnlinked:
#'   To form the \emph{unlinked} crime pairs, \emph{crime groups} are identified 
#'   as the maximal connected offender subgraphs. Then \code{m} indices are drawn
#'   from each crime group (with replacment) and paired with crimes from other crime groups according 
#'   to weights that ensure that large groups don't give the most events.
#'  @examples  
#'  data(crimes)
#'  data(offenders)
#'  seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)
#'  allPairs = makePairs(seriesData,thres=365,m=40)
#'  @export
#'  @name makePairs
##  igraph package is needed for makeUnlinked()
##==============================================================================
NULL


##  makePairs
##==============================================================================
##  Generates indices of linked and unlinked crime pairs (with weights)
##  @return (n x 4) matrix of indices of crime pairs with weights. The last
##    column \code{type} indicates if the crime pair is linked or unlinked.
##  igraph package is needed for makeUnlinked()
#' @rdname makePairs
##==============================================================================
makePairs <- function(X,thres=365,m=40,show.pb=FALSE,seed=NULL){
  linkedPairs = makeLinked(X,thres=thres)          # Get all linked pairs (with weights)
  unlinkedPairs = makeUnlinked(X,m=40,thres=thres,show.pb=show.pb,seed=seed) # Sample unlinked pairs
  linkedPairs$type = 'linked'
  unlinkedPairs$type = 'unlinked'
  allPairs = rbind(linkedPairs,unlinkedPairs)
return(allPairs)
}


##  makeLinked
##==============================================================================
##  Generates unique indices for linked crime pairs (with weights)
##  @return (n x 3) matrix of indices of linked crime pairs with weights.
##  @details The weights are such that each crime series only gives a total weight
##    of 1. Due to co-offending, the sum of weights will be smaller than the
##    number of series with at least two crimes.
##  changed indexing from Index to crimeID columns
#' @rdname makePairs
##==============================================================================
makeLinked <- function(X,thres=365){
  pairwise <- function(A) if(length(A)>1) t(combn(A,2))
  W = tapply(1:nrow(X),X$offenderID,pairwise)
  LP = do.call(rbind,W)  # all linked pairs (indices of X)
  i1 = LP[,1]; i2 = LP[,2]
#  D = data.frame(i1=X$Index[i1],i2=X$Index[i2],offenderID=X$offenderID[i1])
  D = data.frame(i1=X$crimeID[i1],i2=X$crimeID[i2],offenderID=X$offenderID[i1],
                 stringsAsFactors=FALSE)
  #-- threshold based on time
  val = dtdiff(X$TIME[i1],X$TIME[i2],units='days') # time between events
  D = D[val<=thres,]        # remove pairs that exceed threshold
  #- Order indices and remove duplicates
#  flip = D[,1] > D[,2]
  flip = apply(D[,1:2],1,function(x) order(x)[1] != 1L)
  D[flip,1:2] = D[flip,2:1]
  D = unique(D)   # ensure there are no duplicated pairs in a crime series
  #-- Weight each pair according to inverse of number of pairs in series
  tab = as.data.frame(table(D$offenderID))
  ind = match(D$offenderID,tab$Var1)
  EL = data.frame(D,wt=1/tab$Freq[ind]) # linked pairs with weights
  #-- Replace duplicated pairs using smallest weight
  unique.names = paste(EL[,1],EL[,2],sep=":")
  wt = tapply(EL$wt,unique.names,min)
  EL2 = data.frame(EL[match(names(wt),unique.names),1:2],weight=as.numeric(wt))
  return(EL2)
}


##  makeUnlinked
##==============================================================================
##  Generates (a sample of) indices of unlinked crime pairs
##
##  This function generates a set of crimeIDs of unlinked crime pairs. It selects
##  Notes:
##  @details First, \emph{crime groups} are identifyed as the maximal connected
##  offender subgraphs. Then \code{m} indices are drawn from each crime group and
##  paired with crimes from other crime groups according to weights to ensure
##  that large groups don't give the most events.
## igraph package is needed for getting crime groups
## changed indexing from Index to crimeID columns
#' @rdname makePairs
##==============================================================================
makeUnlinked <- function(X,m,thres=365,show.pb=FALSE,seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  #require(igraph) for makeGroups() function
  X$CG = makeGroups(X,method=1)   # get crime groups
  nCG = length(unique(X$CG))      # number of crime groups
  
  Y = tapply(1:nrow(X),X$CG, function(i){
     i = i[!duplicated(X$crimeID[i])] # remove rows with duplicated crimeID
     data.frame(X[i,c('crimeID','CG','TIME'),drop=FALSE],
               wt= 1/((nCG-1)*(length(i))))
  })
  Y = do.call(rbind,Y)
  #-- Sample unlinked events
  getSample <- function(x,...) x[sample.int(length(x), ...)]
  EL = NULL
  I = 1:nrow(Y)
  if(show.pb) pb = txtProgressBar(style=3,max=nCG)
  for(i in 1:nCG){
    ind = which(Y$CG == i)
    i1 = getSample(ind,m,replace=TRUE) # compare m crimes from group i
    i2 = sample(I[-ind],m,prob=Y$wt[-ind]) # sample indices from other groups
    val = dtdiff(Y$TIME[i1],Y$TIME[i2],units='days')
    el = data.frame(i1=Y$crimeID[i1],i2=Y$crimeID[i2],val,stringsAsFactors=FALSE)
    EL = rbind(EL,el)
    if(show.pb) setTxtProgressBar(pb,i)
  }
  if(show.pb) close(pb)
  EL = subset(EL,val<=thres) # EL[,1:2]
#  flip = EL[,1] > EL[,2]
  flip = apply(EL[,1:2],1,function(x) order(x)[1] != 1L)
  EL[flip,1:2] = EL[flip,2:1]
  EL = unique(EL[,1:2])   # unlinked pairs
  EL2 = cbind(EL,weight=1)
  attr(EL2,"num.groups") = nCG   # number of crime groups
return(EL2)
}





##  dtdiff
##==============================================================================
#'  Calculates time between two vectors of datetimes
#'  
##  Inputs:
#'  @param t1 first set of times
#'  @param t2 second set of times
#'  @param units one of: "auto", "secs", "mins", "hours","days", "weeks"
##  Outputs:
#'  @return numeric vector of times between the datetime objects
#'  @keywords internal
##==============================================================================
dtdiff <- function(t1,t2,units='days'){
  as.numeric(abs(difftime(t1,t2,units=units)))
}


## getROC
##==============================================================================
#'  Cacluate ROC like metrics.
#'
#'  Orders scores from largest to smallest and evaluates performance for each
#'    value. This assumes an analyst will order the predicted scores and start
#'    investigating the linkage claim in this order.
##  Inputs:
#'  @param f predicted score for linkage
#'  @param y truth; linked=1, unlinked=0
##  Outputs:
#'  @return data.frame of evaluation metrics:
#'  \itemize{
#'      \item FPR - false positive rate - proportion of unlinked pairs that are
#'        incorrectly assessed as linked
#'      \item TPR - true positive rate; recall; hit rate - proportion of all linked
#'        pairs that are correctly assessed as linked
#'      \item PPV - positive predictive value; precision - proportion of all pairs
#'        that are predicted linked and truely are linked
#'      \item Total - the number of cases predicted to be linked
#'      \item TotalRate - the proportion of cases predicted to be linked
#'      \item threshold - the score threshold that produces the results
#'  }
#'  @examples
#' f = 1:10
#' y = rep(0:1,length=10)
#' getROC(f,y)
#'  @export
##==============================================================================
getROC <- function(f,y){
  ord = order(f,decreasing=TRUE)
  f = f[ord]
  y = y[ord]
  uniq = !duplicated(f,fromLast=TRUE)  # only return results for one duplicated value
  n.pos = sum(y==1)
  n.neg = sum(y==0)
  TP = cumsum(y==1)[uniq]
  FP = cumsum(y==0)[uniq]
  Total = (1:length(f))[uniq]  # total number of cases searched
  TPR = TP/n.pos     # recall
  FPR = FP/n.neg
  PPV = TP/Total     # precision
  data.frame(FPR=FPR,TPR=TPR,PPV=PPV,Total=Total,TotalRate=Total/length(f),
             threshold=f[uniq])
}  
  



