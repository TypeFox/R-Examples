#Get chromosome wise optimal leaf order
#' Chromosome wise leaf ordering
#' 
#' Calculates the optimal leaf ordering pairwise for all linkage groups. 
#'  
#' @param rf Matrix of pairwise recombination frequencies.
#' @param df Vector of cluster numbers, created by splitChr(). Zeros indicated filtered markers and will be ignored.
#' @param method Name of method. Default: seriation (uses the optimal leaf ordering algorithm from the seriation package). 
#' Alternatives endlink (order.endlink from gclus) and endlink-global (ignores linkage groups).
#' @param maxSarf Maximum number of neighbor to include into SARF extension.
#' @return Vector of global marker order.
#' @import stats
#' @examples
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' sortLeafs(rfMat, split)
#' @export
sortLeafs <- function(rf, df, method = "seriation", maxSarf = NULL){
  method <- match.arg(method, c("seriation", "endlink", "endlink-global"))
  if(method == "endlink-global"){
    if (!requireNamespace("gclus", quietly = TRUE)) {
      stop("gclus needed for this function to work. Please install it.",
           call. = FALSE)
    }
    df$order <- gclus::order.endlink(-rf)    
  }else{
    df$order <- rep(0, nrow(rf)) 
    global <- 0
    for(i in 1:max(df$split)){
      sub <- which(df$split == i)
      if(method == "seriation"){
        if (!requireNamespace("seriation", quietly = TRUE)) {
          stop("seriation needed for this function to work. Please install it.",
               call. = FALSE)
        }
#         if (!requireNamespace("gtools", quietly = TRUE)) {
#           stop("gtools needed for this function to work. Please install it.",
#                call. = FALSE)
#         }
        dist <- as.dist(rf[sub, sub])
        hc <- hclust(dist, method = "single")
        # multiple clusters are possible for the same distance, if distances are unique
        hclist <- allTrees(hc, dist, start = 1)
        mode(dist) <- "numeric"
        #seriation::seriate(tmp2, method= "OLO", control = list(hclust=tmp))[[1]]$order
        
        orders <- lapply(hclist, 
                         function(x){
                           hc_obj <-list(merge = x,
                                         method = "single",
                                         order = 1:attr(dist,"Size"))
                           class(hc_obj) <- "hclust"
                           seriation::seriate(dist, method = "OLO", 
                                              control = list(hclust = hc_obj))[[1]]$order
                         } 
                         )
        #orders <- lapply(hclist, function(x) seriation:::.seriate_optimal(list(merge = x), dist)$order)
        #orders <- lapply(hclist, function(x) .Call("order_optimal", dist, x, PACKAGE = "seriation")[[2]])
        if(is.null(maxSarf)){
          maxSarf <- length(sub) - 1
        }
        orders <- do.call(rbind, orders)
        orders <- unique(orders)
        serOut <- sarfExt(orders, dis = dist, maxSarf = maxSarf)
      }else if(method == "endlink"){
        if (!requireNamespace("gclus", quietly = TRUE)) {
          stop("gclus needed for this function to work. Please install it.",
               call. = FALSE)
        }
        serOut <- gclus::order.endlink(-rf[sub, sub])
      }
      df$order[df$split > 0][(global + 1):(global + length(sub))] <- sub[serOut]
      global <- global + length(sub)
    }
  }
  return(df)
}

#' Calculates the SARF value of given input.
#' 
#' The sum of adjecent recombination frequency (SARF) is a measure of how well the marker order is.
#' This function calculates it for a given matrix of pairwise recombination frequencies and marker order.
#' The SARF criterion can be extended to a neighborhood > 1.
#'  
#' @param rf Matrix of pairwise recombination frequencies.
#' @param ord Vector with marker order.
#' @param n Number of neighbors, which are included in the calculation.
#' @return Single numeric value, which is the result of the SARF calculation.
#' @examples
#' data(simTetra)
#' simTetrageno <- bases2genotypes(simTetra, 4)
#' rfMat <- calcRec(simTetrageno, 4)
#' split <- splitChr(rfMat, nchr = 7)
#' split <- sortLeafs(rfMat, split)
#' calcSarf(rfMat, split$order, n = 1)
#' calcSarf(rfMat, split$order, n = 2)
#' calcSarf(rfMat, split$order, n = 3)
#' @references Liu, B.H. 1998, \emph{Statistical genomics: linkage, mapping, and QTL analysis.}
#' @export

calcSarf <- function(rf, ord = 1:(ncol(rf)), n = 1){
  ord <- ord[ord > 0]
  l <- length(ord)
  if(n >= l){
    n <- l - 1
    warning(paste("n was set to ", n))
  }
  rf <- rf[ord, ord]
  sarf <- 0
  for(j in 1:n){
    for(i in (1 + j):l){    
      sarf <- sarf + rf[i, i - j]
    }    
  }
  return(sarf)
}

