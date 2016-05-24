###############################################################
#                                                             #
#   (c) Victor Maus <vwmaus1@gmail.com>                       #
#       Institute for Geoinformatics (IFGI)                   #
#       University of Muenster (WWU), Germany                 #
#                                                             #
#       Earth System Science Center (CCST)                    #
#       National Institute for Space Research (INPE), Brazil  #
#                                                             #
#                                                             #
#   R Package dtwSat - 2015-09-01                             #
#                                                             #
###############################################################


###############################################################
#### TWDTW ALIGNMENT


#' @title Multidimensional Time-Weighted DTW Alignment
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs a multidimensional Time-Weighted DTW 
#' analysis and retrieves the alignments of a query within a time series.
#' 
#' @param query A \link[zoo]{zoo} object with a query time series.
#' @param template A \link[zoo]{zoo} object with a template time series similar 
#' to \code{query}. The \code{template} must have the same number of attributes
#' and be equal to or longer than the \code{query}, 
#' \emph{i.e.} \code{nrow(query)<=nrow(template)}.
#' @param weight A character. ''linear'' for linear weight or ''logistic'' 
#' for logistic weight. Default is NULL that runs the original dtw method,
#' \emph{i.e.} without time weight.
#' @param dist.method A character. Method to derive the local cost matrix.
#' Default is ''Euclidean'' See \code{\link[proxy]{dist}} in package 
#' \pkg{proxy}.
#' @param theta A number. Parameter for ''linear'' time weight. For \code{theta=1} 
#' the time weight is equal to the number of elapsed days.
#' @param alpha A number. The steepness of logistic weight.
#' @param beta A number. The midpoint of logistic weight.
#' @param n.alignments An integer. The maximun number of alignments to 
#' perform. NULL will return all possible alignments. 
#' @param span Span between two points of minimum in days, \emph{i.e.} the minimum  
#' interval between two alignments, for details see [1] 
#' @param step.matrix see \code{\link[dtw]{stepPattern}} in package \pkg{dtw} [2]
#' @param keep preserves the cost matrix, inputs, and other internal structures. 
#' Default is FALSE

#' 
#' @param query.name A query identification.
#' @docType methods
#' @return An object of class \code{\link[dtwSat]{dtwSat-class}}
#' 
#' @references 
#' [1] M\"uller, M. (2007). Dynamic Time Warping. In Information Retrieval for Music 
#' and Motion (pp. 79-84). London: Springer London, Limited. 
#' @references 
#' [2] Giorgino, T. (2009). Computing and Visualizing Dynamic Time Warping Alignments in R: 
#' The dtw Package. Journal of Statistical Software, 31, 1-24.
#' 
#' @seealso \code{\link[dtwSat]{mtwdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' 
#' @examples
#' names(query.list)
#' query.name = "Soybean"
#' alig = twdtw(query.list[[query.name]], template, weight = "logistic", 
#'        alpha = 0.1, beta = 100, n.alignments=4, query.name=query.name)
#' alig
#' 
#' @export
twdtw =  function(query, template, weight=NULL, dist.method="Euclidean",
                  theta=NULL, alpha=NULL, beta=NULL, n.alignments=NULL, span=NULL, 
                  step.matrix = symmetric1, keep=FALSE, query.name=NULL)
{

  if(!is.zoo(query))
    stop("query should be of class zoo.")
  if(!is.zoo(template))
    stop("template should be of class zoo")
  if(ncol(query)!=ncol(template))
    stop("Number of columns in query and in template don't match.")
  if(!is(index(query),"Date"))
    stop("Index in query should be of class Date.")
  if(!is(index(template),"Date"))
    stop("Index in template should be of class Date.")

  .twdtw(query, template, weight, dist.method, theta, alpha, 
         beta, n.alignments, span, step.matrix, keep, query.name)
  
}

#' @title Performs multiple Time-Weighted DTW 
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description This function performs the Time-Weighted DTW for a list 
#' of queries
#' 
#' @param query A \link[zoo]{zoo} object with a query time series
#' @param template A \link[zoo]{zoo} object with a template time series similar 
#' to \code{query}. The \code{template} must have the same number of attributes
#' and be equal to or longer than the \code{query}
#' @param ... additional arguments passed to \code{\link[dtwSat]{twdtw}}
#' @docType methods
#' @return An object of class \link[dtwSat]{dtwSat-class} without \code{internals} 
#' or \code{mapping}
#'
#' @seealso \code{\link[dtwSat]{twdtw}}, \code{\link[dtwSat]{dtwSat-class}}
#' @examples
#' alig = mtwdtw(query.list, template, weight = "logistic", 
#'        alpha = 0.1, beta = 100)
#' alig
#' 
#' @export
mtwdtw = function(query, template, ...){
  if(!is.list(query))
    stop("query is not a list.")
  if(!any(unlist(lapply(query, is, "zoo"))))
    stop("query is not a list of zoo objects.")
  
  query_names = names(query)
  if(is.null(query_names))
    query_names = seq_along(query)
  
  res = new("dtwSat")
  res@internals$template = template
  res@alignments = do.call("rbind", lapply(query_names, function(i){
    alig = twdtw(query=query[[i]], template=template, query.name=i, ...)
    getAlignments(alig)
  }))
  res
}


.twdtw =  function(query, template, weight=NULL, dist.method="Euclidean",
                  theta=NULL, alpha=NULL, beta=NULL, n.alignments=NULL, span=NULL,  
                  step.matrix = symmetric1, keep=FALSE, query.name=NULL)
{

  # Local cost
  phi = proxy::dist(query, template, method=dist.method)
  # Elapsed time
  psi = matrix(0, nrow = nrow(phi), ncol = ncol(phi))
  if(!is.null(weight)){
    psi = .timeCostMatrix(query, template, dist.method)
    psi = switch(weight, 
                 linear   = .linearweight(psi, theta),
                 logistic = .logisticweight(psi, alpha, beta)
    )
  }
  delta = phi + psi
  internals = .costMatrix(delta, step.matrix)
  internals$timeWeight = matrix(psi, nrow = nrow(psi))
  internals$localMatrix = matrix(delta, nrow = nrow(delta))
  internals$query = query
  internals$template = template
  
  alignments = list(quey=numeric(0),from=numeric(0), to=numeric(0), distance=numeric(0), normalizedDistance=numeric(0))
  mapping = list(index1 = numeric(0), index2 = numeric(0))
  d = internals$costMatrix[internals$N,1:internals$M]
  endPoints = .findMin(d, index(template), span=span)
  if(length(endPoints)>0){
    endPoints = endPoints[order(d[endPoints])]
    if(is.null(n.alignments))
      n.alignments = length(endPoints)
    if(length(endPoints) > n.alignments)
      endPoints = endPoints[1:n.alignments]
    # Map low cost paths (k-th paths)
    mapping = lapply(endPoints, function(b){
      return(.kthbacktrack(internals, b))
    })
    # Get the starting point of each path
    startPoints = unlist(lapply(mapping, function(map){
      return(map$index2[1])
    }))

    if(is.null(query.name))
      query.name = 1
    
    alignments = list(query = query.name,
                      from  = index(template)[startPoints],
                      to    = index(template)[endPoints],
                      distance           = d[endPoints],
                      normalizedDistance = d[endPoints] / length(query),                      
                      stringsAsFactors = FALSE)
  }
  
  if(keep) return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping, internals=internals))
  
  return(new("dtwSat", call=match.call(), alignments=alignments, mapping=mapping))
}

.findMin = function(x, timeline, span=NULL){
  NonNA = which(!is.na(x))
  dx = diff(x[NonNA])
  index_min = NonNA[which(dx[-length(dx)] < 0 & dx[-1] >= 0)] + 1
  if(tail(dx,1) < 0)
    index_min = c(index_min,length(x))
  order_min = index_min[order(x[index_min])]
  min_out = array()
  for(i in seq_along(index_min)){
    min_out[i] = order_min[i]
    lower_bound = timeline[order_min[i]] - span
    upper_bound = timeline[order_min[i]] + span
    in_span = lower_bound < timeline[order_min] & timeline[order_min] < upper_bound
    order_min[in_span] = NA
  }
  res = min_out[!is.na(min_out)]
  res
}

.timeCostMatrix = function(query, template, dist.method){ 
  tx = as.numeric(format(index(query), "%j"))
  ty = as.numeric(format(index(template), "%j"))
  psi = proxy::dist(tx, ty, method=dist.method)
  psi[psi>(366/2)] = abs(366 - psi[psi>(366/2)])
  return(psi)
}

.logisticweight = function(x, alpha, beta){
  return( 1 / (1 + exp(1)^(-alpha*(x-beta))) )
}

.linearweight = function(x, theta){
  return( theta * x / 366 )
}

.kthbacktrack = function(alignment, jmin=NULL){
  
  dir = alignment$stepPattern
  npat = attr(dir,"npat")
  
  n = nrow(alignment$costMatrix)
  m = ncol(alignment$costMatrix)
  
  i = n
  j = jmin
  if(is.null(jmin))
    j = alignment$jmin
  
  nullrows = dir[,2]==0 & dir[,3]==0
  tmp = dir[!nullrows,,drop=FALSE]
  
  stepsCache = list()  
  for(k in 1:npat) {
    sbs = tmp[,1]==k  
    spl = tmp[sbs,-1,drop=FALSE]
    nr = nrow(spl)
    stepsCache[[k]] = spl[nr:1,,drop=FALSE]
  }
  
  I = c(i)
  J = c(j)
  
  repeat {
    if(i==1)
      break	
    s = alignment$directionMatrix[i,j]
    if(is.na(s))
      break
    
    steps = stepsCache[[s]]
    ns = nrow(steps)
    
    for(k in 1:ns) {
      if(i-steps[k,1] > 0) {
        I = c(i-steps[k,1],I)
        J = c(j-steps[k,2],J)
      }                         
    }
    
    i = I[1]
    j = J[1]
  }
  
  res = list()
  res$index1 = I
  res$index2 = J
  res
}


.costMatrix = function(delta, step.matrix){
  
  if (!is(step.matrix, "stepPattern"))
    stop("step.matrix is no stepPattern object")
  
  delta = rbind(0, delta)
  n = nrow(delta)
  m = ncol(delta)
  cm = matrix(NA, nrow=n, ncol=m)
  cm[1,] = 0
  
  nsteps = dim(step.matrix)[1]
  cm[1,1] = delta[1,1]
  sm = matrix(NA, nrow=n, ncol=m)

  dir = step.matrix
  npats = attr(dir,"npat")
  for (j in 1:m) {
    for (i in 1:n) {
      if(!is.na(cm[i,j]))
        next
      
      clist = numeric(npats)+NA
      for (s in 1:nsteps) {

        p = dir[s,1]
        I = i-dir[s,2]
        J = j-dir[s,3]
        if(I>=1 && J>=1) {
          cc = dir[s,4]
          if(cc == -1) {
            clist[p] = cm[I,J]
          }else{
            clist[p] = clist[p]+cc*delta[I,J]
          }
        }
      }
      
      minc<-which.min(clist)
      if(length(minc) > 0){
        cm[i,j] = clist[minc]
        sm[i,j] = minc
      }
    }
  }
  
  res = list()
  res$costMatrix = cm[-1,]
  res$directionMatrix = sm
  res$stepPattern = step.matrix
  res$N = n-1
  res$M = m
  res
}







