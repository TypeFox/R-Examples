

# =====================================================
# Function that obtains a statistic of centrality of a
# variable, given a sample of values.
# If the variable is numeric it returns de median, if it
# is a factor it returns the mode. In other cases it
# tries to convert to a factor and then returns the mode.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
centralValue <- function(x, ws=NULL) {
   if (is.numeric(x)) {
     if (is.null(ws)) median(x,na.rm=T)
     else if ((s <- sum(ws)) > 0) sum(x*(ws/s)) else NA
   } else {
     x <- as.factor(x)
     if (is.null(ws)) levels(x)[which.max(table(x))]
     else levels(x)[which.max(aggregate(ws,list(x),sum)[,2])]
   }
 }


# =====================================================
# Small utility function to obtain the target (response)
# variable values corresponding to a given formula and
# data set.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
resp <- function(formula,data) model.response(model.frame(formula, data))


# =====================================================
# Small utility function to obtain the number of the rows
# in a data frame that have either a "large" number of 
# unknown values.
# "Large" can be defined either as a proportion of the
# number of columns or as the number in itself.
# =====================================================
# Luis Torgo, Mar 2009, Mar 2011
# =====================================================
manyNAs <- function(data,nORp=0.2) {
  n <- if (nORp < 1) as.integer(nORp*ncol(data)) else round(nORp,0)
  idxs <- which(apply(data,1,function(x) sum(is.na(x))) > n)
  if (!length(idxs)) warning('Empty index generated, no rows with many NAs. Undesirable effects may be caused if indexing a data frame with this.')
  idxs
}



# =====================================================
# Function that fills in all unknowns using the statistic
# of centrality of the respective column.
# This statistic is either the median for numeric columns
# or the mode for nominal variables.
# =====================================================
# Luis Torgo, Mar 2009
# =====================================================
centralImputation <- function(data) {
  for(i in seq(ncol(data))) 
    if (any(idx <- is.na(data[,i])))
      data[idx,i] <- centralValue(data[,i])
  data
}


# =====================================================
# Function that fills in all unknowns using the k Nearest
# Neighbours of each case with unknows. 
# By default it uses the values of the neighbours and 
# obtains an weighted (by the distance to the case) average
# of their values to fill in the unknows.
# If meth='median' it uses the median/most frequent value,
# instead.
# =====================================================
# Luis Torgo, Mar 2009, Nov 2011
# =====================================================
knnImputation <- function(data,k=10,scale=T,meth='weighAvg',distData=NULL) {

  n <- nrow(data)  
  if (!is.null(distData)) {
    distInit <- n+1
    data <- rbind(data,distData)
  } else distInit <- 1
  N <- nrow(data)

  ncol <- ncol(data)
  nomAttrs <- rep(F,ncol)
  for(i in seq(ncol)) nomAttrs[i] <- is.factor(data[,i])
  nomAttrs <- which(nomAttrs)
  hasNom <- length(nomAttrs)
  contAttrs <- setdiff(seq(ncol),nomAttrs)

  dm <- data
  if (scale) dm[,contAttrs] <- scale(dm[,contAttrs])
  if (hasNom)
    for(i in nomAttrs) dm[,i] <- as.integer(dm[,i])

  dm <- as.matrix(dm)

  nas <- which(!complete.cases(dm))
  if (!is.null(distData)) tgt.nas <- nas[nas <= n]
  else tgt.nas <- nas

  if (length(tgt.nas) == 0)
    warning("No case has missing values. Stopping as there is nothing to do.")

  xcomplete <- dm[setdiff(distInit:N,nas),]
  if (nrow(xcomplete) < k)
    stop("Not sufficient complete cases for computing neighbors.")
  
  for (i in tgt.nas) {

    tgtAs <- which(is.na(dm[i,]))
    
    dist <- scale(xcomplete,dm[i,],FALSE)
    
    xnom <- setdiff(nomAttrs,tgtAs)
    if (length(xnom)) dist[,xnom] <-ifelse(dist[,xnom]>0,1,dist[,xnom])

    dist <- dist[,-tgtAs]
    dist <- sqrt(drop(dist^2 %*% rep(1,ncol(dist))))
    ks <- order(dist)[seq(k)]
    for(j in tgtAs)
      if (meth == 'median')
        data[i,j] <- centralValue(data[setdiff(distInit:N,nas),j][ks])
      else 
        data[i,j] <- centralValue(data[setdiff(distInit:N,nas),j][ks],exp(-dist[ks]))
  }

  data[1:n,]
}




# =====================================================
# Function that inverts the effect of the scale function
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
unscale <- function(vals,norm.data,col.ids) {
  cols <- if (missing(col.ids)) 1:NCOL(vals) else col.ids
  if (length(cols) != NCOL(vals)) stop('Incorrect dimension of data to unscale.')
  centers <- attr(norm.data,'scaled:center')[cols]
  scales <- attr(norm.data,'scaled:scale')[cols]
  unvals <- scale(vals,center=(-centers/scales),scale=1/scales)
  attr(unvals,'scaled:center') <- attr(unvals,'scaled:scale') <- NULL
  unvals
}




# =====================================================
# This function produces a smoothed precision/recall curve
# without the saw-tooth effect of the original curves
# produced by the package ROCR. It is based on the concept
# of interpolated precision
# =====================================================
# Luis Torgo, Feb 2010
# =====================================================
PRcurve <- function(preds,trues,...) {
#  require(ROCR,quietly=T)
  pd <- prediction(preds,trues)
  pf <- performance(pd,'prec','rec')
  pf@y.values <- lapply(pf@y.values,function(x) rev(cummax(rev(x))))
  .plot.performance(pf,...)
}


# =====================================================
# This function produces a cumulative recall chart
# =====================================================
# Luis Torgo, Feb 2010
# =====================================================
CRchart <- function(preds,trues,...) {
#  require(ROCR,quietly=T)
  pd <- prediction(preds,trues)
  pf <- performance(pd,'rec','rpp')
  .plot.performance(pf,...)
}  



# ======================================================================
# Function for normalizing the range of values of a continuous variable.
# Taken from the book "Data preparation for data mining" by Dorian Pyle
# (pp. 271-274)
#
# This function ensures all values will be between 0 and 1.
#
# 13/05/2002, Luis Torgo.
# ----------------------------------------------------------------------
# Example :
# SoftMax(algae[,'NO3'])
# the following obtains the transformation just for one value
# SoftMax(45.23,avg=mean(algae[,'NO3'],na.rm=T),std=sd(algae[,'NO3'],na.rm=T))
#
# Note:
# The lambda parameter controls the range of values that gets a linear
# mapping. It represents the number of standard deviations that should be
# included in the linear mapping region (e.g. 1-> 68% of the distribution gets
# linear mapping, while 2-> 95.5%, 3 -> 99.7%, etc.)
SoftMax <- function(x,lambda=2,avg=mean(x,na.rm=T),std=sd(x,na.rm=T))
{
  if (is.data.frame(x) | is.array(x)) return(apply(x,2,SoftMax,lambda))
  vt <- (x-avg)/(lambda*(std/(2*pi)))
  1/(1+exp(-vt))
}



# ======================================================================
# Function for performing a linear scaling transformation
#
# This function ensures values between 0 to 1 (except for out of the sample
# values, for that use SoftMax).
#
# 13/05/2002, Luis Torgo.
# ----------------------------------------------------------------------
# Example :
# LinearScaling(algae[,'NO3'])
# the following obtains the transformation just for one value
# LinearScaling(45.23,mx=max(algae[,'NO3'],na.rm=T),mn=min(algae[,'NO3'],na.rm=T))
#
LinearScaling <- function(x,mx=max(x,na.rm=T),mn=min(x,na.rm=T))
  (x-mn)/(mx-mn)


# ======================================================================
# Function performs a change of scale
#
# This function ensures values between a given scale
#
# 2/04/2003, Luis Torgo.
# ----------------------------------------------------------------------
#
ReScaling <- function(x,t.mn,t.mx,d.mn=min(x,na.rm=T),d.mx=max(x,na.rm=T)) {
  sc <- (t.mx-t.mn)/(d.mx-d.mn)
  sc*x + t.mn - sc*d.mn
}

