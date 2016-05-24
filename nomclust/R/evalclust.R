#' Evaluation of the Clustering
#' 
#' @description The function evaluates clustering results no matter which clustering method they were obtained by.
#' The clusters are evaluated from a point of view of the within-cluster variability by the following indices:
#' Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
#' Pseudo tau coefficient (PSTau), Pseudo uncertainty coefficient (PSU) and Pseudo F, Indices based on the mutability (PSFM) and the entropy (PSFE).
#' 
#' @param data data frame or matrix with cases in rows and variables in colums. First \code{m1} variables are the
#' original data used for clustering, the next \code{m2} variables express the cluster memberships in an increasing
#' way (e.g. from clu_2 to clu_6).
#' 
#' @param num_var numeric value which determines how many variables in a dataset were used for the clustering.
#' 
#' @param clu_low numeric value expressing the lower bound for number of cluster solutions.
#' 
#' @param clu_high  numeric value expressing the higher bound for number of cluster solutions.
#' 
#' @return Function returns a data frame, where the rows express a serie of cluster solutions and columns
#' clustering evaluation statistics in a following order: \code{WCM}, \code{WCE}, \code{PSTau}, \code{PSU}, \code{PSFM}, \code{PSFE}.
#' 
#' @seealso
#' \code{\link[nomclust]{nomclust}}.
#' 
#' @examples
#' #sample data
#' data(data20)
#' #creation of a dataset with cluster memberships
#' data_clu <- nomclust(data20, iof, clu_high = 7)
#' #binding an original dataset to cluster memberships variables
#' data_clu2 <- cbind(data20, data_clu$mem)
#' #evaluation of created clusters
#' evaluation <- evalclust(data_clu2, 5, clu_high = 7)
#' 
#' @export 

evalclust <- function(data, num_var, clu_low = 2, clu_high = 6) {
  
  if (clu_low >= clu_high) {
    stop("clu_low must be set lower than clu_high")
  }
  
  #if matrix, coerce to data.frame
  data <- as.data.frame(data)
  
  #check lenght of data
  if (ncol(data) != (num_var + length(seq(clu_low:clu_high)))) {
       stop("sum of all set parameters does not match with the size of data")
  }
  
  #cluster membership
  data_clu <- data
  for (i in clu_low:clu_high) {
    names(data_clu)[num_var - clu_low + i + 1] <- paste("clu_", i, sep = "" )
  }
  data <- data_clu[,1:num_var]
  clusters <- data_clu[,(num_var+1):ncol(data_clu)]
  
  #max number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  max_num_cat <- max(num_cat)
  
  #creation of set of 3D matrices
  M <- list()
  for (i in clu_low:clu_high) {
    A <- list()
    A1 <- list()
    MMM <- array(0,dim=c(max_num_cat,i,num_var))
    M1 <- array(0,dim=c(max_num_cat,1,num_var))
    
    for (j in 1:num_var) {
      A[[j]] <- table(data[, j], clusters[,i - clu_low + 1])
      A1[[j]] <- rowSums(A[[j]])
    }
    
    for (j in 1:num_var) {
      MMM[1:nrow(A[[j]]), 1:ncol(A[[j]]), j] <- A[[j]]
      M1[1:nrow(A[[j]]),,j] <- A1[[j]]
    }
    M[[i-clu_low+2]] <- MMM
  }

  #evaluation results
  results <- data.frame(cluster = numeric(clu_high - clu_low + 2), WCM = numeric(clu_high - clu_low + 2), WCE = numeric(clu_high - clu_low + 2),
                        PSTau = numeric(clu_high - clu_low + 2), PSU = numeric(clu_high - clu_low + 2), PSFM = numeric(clu_high - clu_low + 2), PSFE = numeric(clu_high - clu_low + 2))
  
  for (i in clu_low:clu_high) {
    results[i-clu_low+2,1] <- i
    results[i-clu_low+2,2] <- WCM(M[[i-clu_low+2]], num_cat)
    results[i-clu_low+2,3] <- WCE(M[[i-clu_low+2]], num_cat)
    results[i-clu_low+2,4] <- pstau(M[[i-clu_low+2]], M1, num_cat)
    results[i-clu_low+2,5] <- psu(M[[i-clu_low+2]], M1, num_cat)
    results[i-clu_low+2,6] <- psfm(M[[i-clu_low+2]], M1, num_cat)
    results[i-clu_low+2,7] <- psfe(M[[i-clu_low+2]], M1, num_cat)
    results[1,1] <- 1
    results[1,2] <- WCM(M1, num_cat)
    results[1,3] <- WCE(M1, num_cat)
    results[1,4:7] <- NA
  }
  return(results)
}