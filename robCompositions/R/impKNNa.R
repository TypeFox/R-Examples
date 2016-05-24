#' Imputation of missing values in compositional data using knn methods
#' 
#' This function offers several k-nearest neighbor methods for the imputation
#' of missing values in compositional data.
#' 
#' The Aitchison \code{metric} should be chosen when dealing with compositional
#' data, the Euclidean \code{metric} otherwise.
#' 
#' If \code{primitive} \eqn{==} FALSE, a sequential search for the
#' \eqn{k}-nearest neighbors is applied for every missing value where all
#' information corresponding to the non-missing cells plus the information in
#' the variable to be imputed plus some additional information is available. If
#' \code{primitive} \eqn{==} TRUE, a search of the \eqn{k}-nearest neighbors
#' among observations is applied where in addition to the variable to be
#' imputed any further cells are non-missing.
#' 
#' If \code{normknn} is TRUE (prefered option) the imputed cells from a nearest
#' neighbor method are adjusted with special adjustment factors (more details
#' can be found online (see the references)).
#' 
#' @param x data frame or matrix
#' @param method method (at the moment, only \dQuote{knn} can be used)
#' @param k number of nearest neighbors chosen for imputation
#' @param metric \dQuote{Aichison} or \dQuote{Euclidean}
#' @param agg \dQuote{median} or \dQuote{mean}, for the aggregation of the
#' nearest neighbors
#' @param primitive if TRUE, a more enhanced search for the $k$-nearest
#' neighbors is obtained (see details)
#' @param normknn An adjustment of the imputed values is performed if TRUE
#' @param das depricated. if TRUE, the definition of the Aitchison distance,
#' based on simple logratios of the compositional part, is used (Aitchison,
#' 2000) to calculate distances between observations.  if FALSE, a version
#' using the clr transformation is used.
#' @param adj either \sQuote{median} (default) or \sQuote{sum} can be chosen
#' for the adjustment of the nearest neighbors, see Hron et al., 2010.
#' @return \item{xOrig }{Original data frame or matrix} \item{xImp }{Imputed
#' data} \item{w }{Amount of imputed values} \item{wind }{Index of the missing
#' values in the data} \item{metric }{Metric used}
#' @author Matthias Templ
#' @export
#' @seealso \code{\link{impCoda}}
#' @references Aitchison, J. and Barcelo-Vidal, C. and Martin-Fernandez, J.A.
#' and Pawlowsky-Glahn, V. (2000) Logratio analysis and compositional distance,
#' Mathematical Geology 32(3):271-275.
#' 
#' Hron, K. and Templ, M. and Filzmoser, P. (2010) Imputation of missing values
#' for compositional data using classical and robust methods
#' \emph{Computational Statistics and Data Analysis}, vol 54 (12), pages
#' 3095-3107.
#' @keywords manip multivariate
#' @examples
#' 
#' data(expenditures)
#' x <- expenditures
#' x[1,3]
#' x[1,3] <- NA
#' xi <- impKNNa(x)$xImp
#' xi[1,3]
#' 
`impKNNa` <-
function(x, method="knn", k=3,
                    metric="Aitchison", agg="median", primitive=FALSE, normknn=TRUE, das=FALSE, adj="median"){
### Coda Version 2!
### Nearest Neighbor imputation algorithm of Missing values
### based on distance matrices.
###
### Program by Matthias Templ, Oktober 2008
### R Version: 2.7.1
### Copyright : TU-Wien, 2008
###

# TODO: help file mit adj ergaenzen


x <- as.matrix(x)
class(x) <- "matrix"

N <- dim(x)[1]
P <- dim(x)[2]

if(any(rowSums(is.na(x)) == P)) stop("One or more observations constist of only NA's")

xOrig <- xmiss <- x

w <- is.na(x)
w2 <- !is.na(x)
if(metric=="Euclidean" & primitive==FALSE){
findknn <- function(x, i, j){
  ## find knn
  m1 <- which(!is.na(x[i,]))
  mi <- which(!is.na(x[,j,drop=FALSE]) & !is.na(rowSums(x[,m1,drop=FALSE])) )
  d <- rowSums(t((t(x[mi,m1,drop=FALSE]) - x[i,m1])^2))
  names(d) <- mi
  as.numeric(names(which(d <= quantile(d, k/N))))
}
}
if(metric=="Euclidean" & primitive==TRUE){
findknn <- function(x, i, j){
  ## find knn
  m1 <- which(!is.na(x[i,]))
  mi <- which(!is.na(x[,j,drop=FALSE]) )
  d <- rowSums(t((t(x[mi,m1,drop=FALSE]) - x[i,m1])^2),na.rm=TRUE)
  names(d) <- mi
  as.numeric(names(which(d <= quantile(d, k/N))))
}
}

if(metric=="Aitchison" & primitive==FALSE){
findknn <- function(x, i, j){
  m1 <- which(!is.na(x[i,]))
  mi <- which(c(!is.na(x[,j,drop=FALSE])) & c(!is.na(rowSums(x[,m1,drop=FALSE]))))
  xclr <- cenLR(rbind(x[mi, m1, drop=FALSE], x[i, m1]))$x.clr

  d <- rowSums(t(abs(t(xclr[-nrow(xclr),]) - c(xclr[nrow(xclr),]))))
  names(d) <- mi
  wA <- which(d <= quantile(d, k/length(d)))
### new change for zeros paper:
#			print(is.numeric(xclr))
#if(is.numeric(xclr)) xclr <- matrix(xclr, ncol=length(xclr))
#d <- rowSums(t(abs(t(xclr[-nrow(xclr),,drop=FALSE]) - c(xclr[nrow(xclr),,drop=FALSE]))), na.rm=TRUE)
#names(d) <- mi
#print(d)
#wA <- which(d <= quantile(d, k/length(d)))
### end new change
  w <- as.numeric(names(wA))
  list(knn=w, m1=m1)
}
}
if(metric=="Aitchison" & primitive==TRUE & das == TRUE){
findknn <- function(x, i, j){
  m1 <-  which(!is.na(x[i,]))
  mi <- which(c(!is.na(x[,j,drop=FALSE])) )
    da <- function(x,y){
    d <- 0
    p <- length(x)
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        d <- d + (log(x[i]/x[j]) - log(y[i]/y[j]))^2
      }
    }
    d=d/p
    d
  }
  ref <- x[i,m1]
  d <- apply(x[mi, m1, drop=FALSE], 1, function(z){da(x=z, y=ref)})
  names(d) <- mi
  wA <- which(d <= quantile(d, k/length(d)))  
  w <- as.numeric(names(wA))
  list(knn=w, m1=m1)      
}
}

if(metric=="Aitchison" & primitive==FALSE & das == TRUE){
findknn <- function(x, i, j){
  m1 <-  which(!is.na(x[i,]))
  mi <- which(c(!is.na(x[,j,drop=FALSE])) & c(!is.na(rowSums(x[,m1,drop=FALSE]))))
    da <- function(x,y){
    d <- 0
    p <- length(x)
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        d <- d + (log(x[i]/x[j]) - log(y[i]/y[j]))^2
      }
    }
    d=d/p
    d
  }
  ref <- x[i,m1]
  d <- apply(x[mi, m1, drop=FALSE], 1, function(z){da(x=z, y=ref)})
  names(d) <- mi
  wA <- which(d <= quantile(d, k/length(d)))  
  w <- as.numeric(names(wA))
  list(knn=w, m1=m1)      
}
}
if(metric=="Aitchison" & primitive==TRUE){
findknn <- function(x, i, j){
  m1 <- which(!is.na(x[i,]))
  a <- c(!is.na(x[,j,drop=FALSE]))
  b <- ifelse(rowSums(is.na(x[,-j])) == 0, TRUE, FALSE)
  mi <- which(a & b)
  xclr <- cenLR(rbind(x[mi, m1, drop=FALSE], x[i, m1]))$x.clr
  d <- rowSums(t(abs(t(xclr[-nrow(xclr),]) - c(xclr[nrow(xclr),]))), na.rm=TRUE)
  #d[d == 0] <- NA   ## dirty
  names(d) <- mi
  wA <- which(d <= quantile(d, k/length(d), na.rm=TRUE))
  w <- as.numeric(names(wA))
  list(knn=w, m1=m1)
}
}

if(metric=="Euclidean"){
imp <- function(x, i, j){
  ## workhorse: do the imputation
  knn <- findknn(x, i, j)
  get(agg)(x[knn,j], na.rm=TRUE) #median
}
}

if(metric=="Aitchison" & normknn == FALSE){
imp <- function(x, i, j){
  ## workhorse: do the imputation
  knn <- findknn(x, i, j)
  ### median, normalization:
  #fac <- knn$xclriSum/knn$medKnnSum
  #median(x[knn$knn, j])#*fac  # fac anders berechnen
  medSum <- get(adj)(apply(x[knn$knn, knn$m1, drop=FALSE], 2, function(x,...){get(agg)(x, na.rm=TRUE)} ), na.rm=TRUE) 
  xSum <- get(adj)(x[i, knn$m1, drop=FALSE], na.rm=TRUE)
  fac <- xSum/medSum
  get(agg)(x[knn$knn, j], na.rm=TRUE)*fac #median
}
}

if(metric=="Aitchison" & normknn == TRUE){
imp <- function(x, i, j){
  knn <- findknn(x, i, j)
  ##normalization:
  #xSum <- sum(x[i, knn$m1, drop=FALSE], na.rm=TRUE)
  xSum <- get(adj)(x[i, knn$m1, drop=FALSE], na.rm=TRUE)
  if( length(knn$knn) == 1){
    #knnSum <- sum(x[knn$knn, knn$m1], na.rm=TRUE)
    knnSum <- get(adj)(x[knn$knn, knn$m1], na.rm=TRUE) 
  } else{
    #knnSum <- rowSums(x[knn$knn, knn$m1, drop=FALSE], na.rm=TRUE)
    knnSum <- apply(x[knn$knn, knn$m1, drop=FALSE], 1, function(x,...){get(agg)(x, na.rm=TRUE)} )
  }
  fac <- xSum/knnSum
  get(agg)(x[knn$knn, j] * fac, na.rm=TRUE)
  ## statt Mean das geom. Mittel?
}
}

#if(metric == "Aitchison" & ratios == TRUE){
#    knn <- findknn(x, i, j)
#    p <- ncol(x)
#    y <- apply(apply(x[knn$knn, , drop = FALSE], 1, function(x){x[2:p]/x[1]}),1,median)
#    x[i,j] <- 1
#    x[i,-j] <- y
#}

#if(metric=="Aitchison" & unweighted == TRUE){
#  imp <- function(x, i, j){
#    knn <- findknn(x, i, j)
#    x[i,j] <- get(agg)(x[knn$knn, j])
#  }  
#}

for(i in 1:N){
  for(j in 1:P){
    if(is.na(x[i,j])){
      x[i,j] <- imp(xmiss, i, j)
    }
  }
}

w <- is.na(xOrig)
colnames(x) <- colnames(xOrig)
res <- list(xOrig=xOrig, xImp=x, criteria=NULL, iter=NULL, w=length(which(w)), wind=w, metric=metric)
class(res) <- "imp"
res

}

