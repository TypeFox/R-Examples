#' Profile Analysis via Multidimensional Scaling
#' 
#' The \code{pams} function implements profile analysis via multidimensional scaling as described by Davison, Davenport, and Bielinski (1995) and Davenport, Ding, and Davison (1995).
#' 
#' The \code{pams} function computes similarity/dissimilarity indices based on Euclidean distances between the scores provided in the data, and then extracts dimensional coordinates for each score using multidimensional scaling. A weight matrix, level parameters, and fit measures are computed for each subject in the data. 
#' 
#' @export
#' @importFrom stats cmdscale dist
#' 
#' @param data A data matrix or data frame; rows represent individuals, columns represent scores; missing scores are not allowed.
#' @param dim Number of dimensions to be extracted from the data.
#' 
#' @return 
#' \itemize{
#' \item \code{dimensional.configuration} - A matrix that provides prototypical profiles of dimensions extracted from the data.
#' \item \code{weights.matrix} - A matrix that includes the subject correspondence weights for all dimensions, level parameters, and the subject fit measure which is the proportion of variance in the subject's actual profiles accounted for by the prototypical profiles.}
#' 
#' @references Davenport, E. C., Ding, S., & Davison, M. L. (1995). PAMS: SAS Template. 
#' @references Davison, M. L., Davenport, E. C., & Bielinski, J. (1995). PAMS: SPSS Template. 
#' @examples 
#' \dontrun{
#' data(PS)
#' result <- pams(PS[,2:4], dim=2)
#' result
#' }
#' @seealso \code{\link{cpa}}, \code{\link{pr}}
#' 

pams <- function(data, dim) {
  
  #PART I
  rawdata <- as.matrix(data)
  raw <- t(rawdata)
  k <- nrow(raw)
  
  distfile <- as.matrix(dist(raw, method = "euclidean",))
  coord <- cmdscale(distfile, k = dim)
  coord <- apply(coord,2,function(x) round(((x/k)*(-1)), digits=6))
  
  
  #PART II
  colnames(coord) <- c(paste("Dimension",1:dim,sep=""))
  m <- as.matrix(rawdata)
  colnames(m) <- c(paste("scale",1:k,sep=""))
  y <- coord
  
  R  <- nrow(y)
  COL <- matrix(1,nrow=R,ncol=1)
  y1 <- cbind(y,COL)
  M1 <- t(y1)%*%(y1)
  M2 <- t(y1)%*%t(m)
  W  <- solve(M1,M2)
  TW <- t(W)
  
  
  #PART III
  m1 <- TW%*%t(y1)
  k <- ncol(m)
  r <- nrow(m)
  col <- matrix(1,nrow=1,ncol=k)
  m1.1 <- as.matrix(rowSums(m1),ncol=1)
  mlrsum <- m1 - ((m1.1%*%col)/k)
  pvar <- apply(mlrsum,2,function(x) x^2)%*%t(col)
  m1.2 <- as.matrix(rowSums(m),ncol=1)
  mrsum <- m - ((m1.2%*%col)/k)
  var <- apply(mrsum,2,function(x) x^2)%*%t(col)
  col <- pvar/var
  
  
  #PART IV
  w1 <- cbind(TW,col)
  weights.matrix1 <- apply(w1[,1:dim],2,function(x) round(x, digits=6))
  weights.matrix2 <- apply(w1[,(dim+1):(dim+2)],2,function(x) round(x, digits=2))
  weights.matrix  <- cbind(weights.matrix1,weights.matrix2)
  colnames(weights.matrix) <- c(paste("weight",1:dim,sep=""),"level","R.sq")
  
  output <- list(weights.matrix=weights.matrix, dimensional.configuration=coord)
  
  return(output)
}

