#' Nominal Clustering
#' 
#' @description The Nominal Clustering (nomclust) performs hierarchical cluster analysis (HCA) with objects characterized by nominal (categorical) variables.
#' It performs a serie of cluster solutions, usually from two-cluster solution till six-cluster solution.
#' It allows to choose one from 11 different similarity measures and one from 3 linkage methods. The function
#' also contains an evaluation part. The created clusters are evaluated from a point of view of the
#' within-cluster variability by the following indices: Within-cluster Mutability coefficient (\code{WCM}),
#' Within-cluster entropy coefficient (\code{WCE}), Pseudo tau coefficient (\code{PSTau}), Pseudo uncertainty coefficient (\code{PSU}) and Pseudo F
#' Indices based on the mutability (\code{PSFM}) and the entropy (\code{PSFE}).
#' 
#' @param data data frame or a matrix with cases in rows and variables in colums. Cases are characterized by nominal (categorical) variables coded as numbers.
#' 
#' @param measure character string defining the similarity measure which wil be used for computation of proximity matrix:
#' \code{"eskin"}, \code{"good1"}, \code{"good2"}, \code{"good3"}, \code{"good4"}, \code{"iof"}, \code{"lin"}, \code{"lin1"}, \code{"morlini"}, \code{"of"}, \code{"sm"}.
#' 
#' @param method character string defining the clustering method. The following methods can be used: \code{"average"}, \code{"complete"}, \code{"single"}.
#' 
#' @param clu_low  numeric value expressing the lower bound for number of cluster solutions.
#' 
#' @param clu_high  numeric value expressing the higher bound for number of cluster solutions.
#' 
#' @param eval logical operator; if TRUE, there is performed an evaluation of clustering results
#' 
#' @param prox logical operator; if TRUE, the proximity matrix is a part of the output
#' 
#' @return Function returns a list following components:
#' \cr
#' \cr
#' \code{mem}  data frame consisting of cluster membership variables 
#' \cr
#' \cr
#' \code{eval} data frame containing clustering evaluation statistics
#' \cr
#' \cr
#' \code{prox} matrix containing proximities between all combination of pairs of objects (voluntary)
#'
#' @seealso
#' \code{\link[nomclust]{evalclust}}, \code{\link[cluster]{agnes}}.
#' 
#' @examples
#' #sample data
#' data(data20)
#' hca <- nomclust(data20, iof, method = "average", clu_high = 5, prox = TRUE)
#' #getting evaluation statistics
#' eval <- hca$eval
#' #getting cluster membership variables
#' mem <- hca$mem
#' #getting a proximity matrix
#' prox <- hca$prox
#' 
#' @export

nomclust <- function (data, measure = iof, clu_low = 2, clu_high = 6, eval = TRUE, 
                      prox = FALSE, method = "complete") {
  
  #packages needed
  #require(cluster)
  #check for cluster sizes
  if (clu_low >= clu_high) {
    stop("clu_low must be set lower than clu_high")
  }
  #elimination of other clustering methods
  if (method %in% c("ward", "weighted", "flexible", "gaverage")) {
    stop("invalid clustering method")
  }
  
  #if matrix, coerce to data.frame
  if(is.matrix(data) == 1) {
    data <- as.data.frame(data)
  }
  
  #computing the proximity matrices
  prox_matrix <- measure(data)
  
  #number of variables of dataset
  num_var <- ncol(data)
  
  #max number of categories
  num_cat <- sapply(data, function(x) length(unique(x)))
  max_num_cat <- max(num_cat)

  #hierarchical cluster analysis, where "prox" is a proximity matrix
  hca <- agnes(prox_matrix, diss = TRUE, method = method)

  #cluster membership
  data_clu <- data
  for (i in clu_low:clu_high) {
   clusters <- cutree(hca, i)
   data_clu <-data.frame(data_clu, clusters)
   names(data_clu)[num_var - clu_low + i + 1] <- paste("clu_", i, sep = "" )
  }
  clusters <- data_clu[,(num_var+1):ncol(data_clu)]
  
  if (eval == 1) {
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
  }

  if (eval == 1 & prox == 1) {
    list <- list(mem = clusters, eval = results, prox = prox_matrix)
  }
  if (eval == 0 & prox == 1) {
    list <- list(mem = clusters, prox = prox_matrix)
  }
  if (eval == 1 & prox == 0) {
    list <- list(mem = clusters, eval = results)
  }
  if (eval == 0 & prox == 0) {
    list <- list(mem = clusters)
  }

  return(list)
}
  