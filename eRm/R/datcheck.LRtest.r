datcheck.LRtest <- function(x, X, model)
{
#sanity checks for LRtest (internal function of LRtest.R)
#x...submatrix (splitted with "splitcr" and called within Xlist) 
#X...original data matrix (from model fit)

exclude <- NULL                                             #vector with items to be excluded

#----check full/0 responses------
n.NA <- colSums(apply(X,2,is.na))                                   #number of NA's per column
maxri <- (dim(X)[1]*(apply(X,2,max,na.rm=TRUE)))-n.NA               #maximum item raw scores with NA
ri <- apply(x,2,sum,na.rm=TRUE)                              #item raw scores
exclude <- c(exclude,which((ri==maxri) | (ri==0)))  

#----check full(-1) NA's---------
allna.vec <- apply(x,2,function(y) {
                         naTF <- is.na(y)
                         (sum(naTF) >= length(y-1)) 
                         }) 
exclude <- c(exclude,which(allna.vec))

#----minimum category = 0--------
ri.min <- apply(x,2,min,na.rm=TRUE)                                 #if no 0 responses
exclude <- c(exclude,which(ri.min!=0))

#----RSM-checks for same number of categories--------
if ((model == "RSM") || (model == "LRSM")) {
   highcat <- max(X, na.rm=TRUE)                    #highest category in original data
   highcat.sub <- apply(x,2,max,na.rm=TRUE)             #RSM check for equal number of categories
   exclude <- c(exclude,which(highcat.sub != highcat))
}

#---PCM checks for all categories responses---------
if ((model=="PCM") || (model=="LPCM")) {                         #check if there are missing categories for PCM (for RSM doesn't matter)
  cat.data <- apply(X,2,function(y) list(unique(na.exclude(y)))) #categories of orginal data
  cat.sub <- apply(x,2,function(y) list(unique(na.exclude(y))))  #categories of subgroup data
  catcomp <- mapply(function(y.s,y.d) {
                      (length(y.s[[1]]) == (length(y.d[[1]])))
                    },cat.sub,cat.data)
  exclude <- c(exclude,which(!catcomp))
}

return(unique(exclude))             #return vector with items to be eliminated
}