#' compute win-loss probabilities
#' 
#' \code{conductance} compute win-loss probabilities for all possible pairs
#'  based upon the combined information from directed wins/losses and 
#'  indirect win/loss pathways from the network.
#' 
#' @param conf a matrix of conf.mat class. An N-by-N conflict matrix whose \code{(i,j)}th element is the number of times i defeated j.
#' @param maxLength an integer greater than 1 and less than 7, indicating the maximum length of paths to identify. 
#' @param alpha a positive integer that 
#' reflects the influence of an observed win/loss interaction 
#' on an underlying win-loss probability. 
#' It is used in the calculation of the posterior distribution 
#' for the win-loss probability of \code{i} over \code{j}: \eqn{Beta(\alpha c_{i,j} +\beta, c_{i,j}+\beta)}. 
#' In the absence of expertise to accurately estimate alpha, 
#' it is estimated from the data.
#' @param beta a positive numeric value that, like alpha, 
#' reflects the influence of an observed win/loss interaction 
#' on an underlying win-loss probability. 
#' Both \eqn{\alpha} and \eqn{\beta} are chosen such that \eqn{((\alpha + \beta)/(\alpha + 2\beta))^2} is 
#' equal to the order-1 transitivity of the observed network. 
#' Therefore, \eqn{\beta} is commonly set to 1.
#' @param strict a logical vector of length 1. It is used in transitivity definition for alpha estimation. 
#' It should be set to TRUE when a transitive triangle is defined as all pathways in the triangle go to the same direction;
#' it should be set to FALSE when a transitive triangle is defined as PRIMARY pathways in the triangle go to the same direction.
#' Strict = FALSE by default.
#' @return a list of two elements. 
#' 
#'  \item{imputed.conf}{An N-by-N conflict matrix whose \code{(i,j)}th element is the 
#'    'effective' number of wins of \code{i} over \code{j}.}
#'    
#'  \item{p.hat}{An N-by-N numeric matrix whose \code{(i,j)}th element is the estimated 
#'      win-loss probability. 
#'      Three functions (\code{\link{valueConverter}}, \code{\link{individualDomProb}}, and \code{\link{dyadicLongConverter}}) are provided to convert win-loss probability 
#'      into other formats that are easier for further analysis of win-loss probability. }
#'      
#' @details This function performs two major steps. 
#' First, repeated random walks through the empirical network 
#' identify all possible directed win-loss pathways 
#' between each pair of nodes in the network. 
#' Second, the information from both direct wins/losses and 
#' pathways of win/loss interactions are combined into an estimate of 
#' the underlying probability of \code{i} over \code{j}, for all \code{ij} pairs.
#' 
#' @seealso \code{\link{as.conflictmat}}, \code{\link{findIDpaths}}, \code{\link{transitivity}}, \code{\link{simRankOrder}}
#' 
#' @references Fushing H, McAssey M, Beisner BA, McCowan B. 2011. 
#' Ranking network of a captive rhesus macaque society: a sophisticated corporative kingdom. 
#' PLoS ONE 6(3):e17817.
#' 
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find win-loss probability matrix
#' perm2 <- conductance(confmatrix, 2, strict = FALSE)
#' perm2$imputed.conf
#' perm2$p.hat
#' @export

conductance = function(conf, maxLength, alpha = NULL, beta = 1, strict = FALSE){
  
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf))){
    conf = as.conflictmat(conf)
  }
  
  if(maxLength < 2 | maxLength > 6) {
    stop("'maxLength' should be an integer greater than 1 and less than 7.")
  }
  
  # making sure maxLength is an integer
  if(maxLength %% as.integer(maxLength) != 0) {
    stop("'maxLength' needs to be an integer.")
  }
  
  N = nrow(conf)
  
  ### percMat will contain direct + indirect information from win-loss paths
  percMat = conf
  
  outdegree = rowSums(conf)
  
  # calculate alpha if not exist
  conf.trans <- Perc::transitivity(conf, strict = FALSE)
  if (is.null(alpha)) {
    alpha <- conf.trans$alpha
  }
  # if alpha is larger than 500, use 500.
  alpha <- min(alpha, 500)
  
  paths = allPaths(conf, maxLength)
  
  ### Populating the direct + indirect conflict matrix called "percMat"
  
  # Aaron's update
# temp disable
#    for(k in 1:(maxLength - 1)){
#    for(r in 1:nrow(paths[[2]][[k]])){
#      wij.star = diag(conf[paths[[2]][[k]][r, 1:(k+1)], 
#                           paths[[2]][[k]][r, 2:(k+2)]])
#      wji.star = diag(conf[paths[[2]][[k]][r, 2:(k+2)], 
#                           paths[[2]][[k]][r, 1:(k+1)]])
#      percMat[paths[[2]][[k]][r,1], paths[[2]][[k]][r,k+2]] = 
#        percMat[paths[[2]][[k]][r,1], paths[[2]][[k]][r,k+2]] + 
#        prod((alpha * wij.star + beta) / 
#               (alpha * (wij.star + wji.star) + 2*beta)/mean(outdegree))
#    }
#  }
  
# fix codes which used pathways incorrectly  
  if(sum(conf[row(conf) != col(conf)] == 0) > 0){
    paths = allPaths(conf, maxLength)
    
    ### Populating the direct + indirect conflict matrix called "percMat"
    for(k in 1:(maxLength - 1)){
      for(r in 1:nrow(paths[[2]][[k]])){
        percMat[paths[[2]][[k]][r,1], paths[[2]][[k]][r,k+2]] = 
          percMat[paths[[2]][[k]][r,1], paths[[2]][[k]][r,k+2]] + 
          ((alpha + beta)/(alpha + 2*beta)/mean(outdegree))^k
        # gc()
      }
    }
  }
  else{
    constant = sum(((alpha + beta)/(alpha + 2 * beta)/mean(outdegree))
                   ^(1:(maxLength - 1)))
    percMat = percMat + constant
    diag(percMat) = 0
  }
  
  
  ### "percMat2" is the estimated win-loss probability matrix
  
  percMat2 = matrix(0, N, N)
  for(i in 2:N){
    for(j in 1:(i-1)){
      temp1 = (alpha * percMat[i,j] + beta)/(alpha * percMat[i,j] +  # percMat[i, j]: times i triumphs j
                                               alpha * percMat[j,i] + 2 * beta) # percMat[j, i]: times j triumphs i.
      temp2 = (alpha * percMat[j,i] + beta)/(alpha * percMat[i,j] + 
                                               alpha * percMat[j,i] + 2 * beta)
      percMat2[i,j] = ifelse(is.nan(temp1), 0.5, temp1)  # lower triangle.
      percMat2[j,i] = ifelse(is.nan(temp2), 0.5, temp2)  # corresponding upper triangle
      # if(verbose){print(i)}  # "Error: object 'verbose' not found"
    }
  }
  row.names(percMat2) <- row.names(conf)
  colnames(percMat2) <- colnames(conf)
  return(list(imputed.conf = percMat, p.hat = percMat2))  
}

#' win-loss probability matrix value converter
#' 
#' \code{valueConverter} converts or transforms all values (which range from 0.0 to 1.0)
#'  in the win-loss probability matrix into 0.5 - 1.0
#' 
#' @param matrix the win-loss matrix which is the second output from \code{conductance}. 
#' @return a matrix of win-loss probability ranging from 0.5 - 1.0.
#' @seealso \code{\link{conductance}}, \code{\link{individualDomProb}}, \code{\link{dyadicLongConverter}}
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find win-loss probability matrix
#' perm2 <- conductance(confmatrix, 2)
#' perm2$imputed.conf
#' perm2$p.hat
#' convertedValue <- valueConverter(perm2$p.hat)
#' @export

valueConverter <- function(matrix){
  
  if (!(is.matrix(matrix))) {
    stop("Only matrix is accepted as input.")
  }
  
  matrixAbove0.5 <- abs(0.5 - matrix) + 0.5
  return(matrixAbove0.5)
}

#' dyadic long format converter
#' 
#' \code{dyadicLongConverter} convert win-loss probability matrix into long format for each dyad
#' 
#' @param matrix the win-loss matrix which is the second output from \code{conductance}. 
#' @return a dataframe of dyadic level win-loss probability and ranking certainty.
#' @details values on the diagonal of the matrix are not included in the converted long-format data.
#' @seealso \code{\link{conductance}}, \code{\link{valueConverter}}, \code{\link{individualDomProb}}
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find win-loss probability matrix
#' perm2 <- conductance(confmatrix, 2)
#' perm2$imputed.conf
#' perm2$p.hat
#' dl <- dyadicLongConverter(perm2$p.hat)
#' @export

# to do:
# 1. change ID2 into character vector; make sure both IDs are character vector.
# 2. clarify in documentation that each unique dyad appear only once.
# 3. clarify in documentation that which ID is the dominant one.
dyadicLongConverter <- function(matrix){
  if (!(is.matrix(matrix))) {
    stop("Only matrix is accepted as input.")
  }
  matrix[lower.tri(matrix, diag = TRUE)] <- NA
  dp.df <- as.data.frame(matrix)
  dp.df2 <- dp.df
  dp.df2$rowID <- rownames(dp.df)
  dp.long <- reshape2::melt(dp.df2, 
                            id.vars = "rowID", 
                            variable.name = "ID2",
                            value.name = "ID1 Win Probability")
  names(dp.long)[1] <- "ID1"
  dpComplete <- dp.long[complete.cases(dp.long), ]
  dpComplete[,"ID2 Win Probability"] <- 1 - dpComplete[,3]
  dpComplete$RankingCertainty <- abs(0.5 - dpComplete[,4]) + 0.5
  return(dpComplete)
}


#' individual-level probability converter
#' 
#' \code{individualDomProb} convert win-loss probability matrix into long format for each dyad
#' 
#' @param matrix the win-loss matrix which is the second output from \code{conductance}. 
#' @return a dataframe. Averaging probability of win-loss relationship with all other individuals.
#' 
#' @seealso \code{\link{conductance}}, \code{\link{valueConverter}}, \code{\link{dyadicLongConverter}}
#' 
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find win-loss probability matrix
#' perm2 <- conductance(confmatrix, 2)
#' perm2$imputed.conf
#' perm2$p.hat
#' individualLevelOutput <- individualDomProb(perm2$p.hat)
#' @export

individualDomProb <- function(matrix){
  
  if (!(is.matrix(matrix))) {
    stop("Only matrix is accepted as input.")
  }
  
  matrixAbove0.5 <- valueConverter(matrix)
  
  Mean <- apply(data.frame(valueConverter(matrixAbove0.5)), 2, mean)
  SD <- apply(data.frame(valueConverter(matrixAbove0.5)), 2, sd)
  
  attributes(Mean) <- NULL
  attributes(SD) <- NULL
  
  individualProb <- data.frame(ID = rownames(matrix), Mean = Mean, SD = SD)
  return(individualProb)
}


# to do:
# -- more explanations for alpha. to add - allow user to set alpha
# -- more explanations for beta