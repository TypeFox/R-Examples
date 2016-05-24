

#' Computes phylogenetic signal with different methods
#'
#' This function computes phylogenetic signal statistics
#' (Blomberg's K and K*, Abouheif's Cmean, Moran's I, and Pagel's Lambda) for traits in \code{phylo4d} objects.
#'
#' @param p4d a \code{phylo4d} object.
#' @param methods a character vector giving the methods to compute phylogenetic signal (see Details).
#' @param reps an integer. The number of repetitions for the estimation of p.values with randomization.
#' @param W an optional matrix of phylogenetic weights to compute Moran's I. By default the matrix
#' is computed with the function \code{\link[adephylo]{proxTips}} with patristic distances.
#'
#'@details p4d must be a \code{phylo4d} object as defined in \pkg{phylobase} package.
#'By default, the \code{methods} argument is set to "\code{all}" and all the available methods are used.
#'The user can specify which method(s) to use. Possible methods are
#'"\code{I}" (Gittleman & Kot 1990), "\code{Cmean}" (Abouheif 1999)", "\code{Lambda}" (Pagel 1999),
#'"\code{K}" and "\code{K.star}" (Blomberg et al. 2003).
#'
#'@return A list of two dataframes with the values of statistics and associated
#'p.values for each tested trait and method.
#'
#' @author This function is a general wrapper for C++ subroutines.
#' C++ code is adapted from R functions in Pavoine and Ricotta (2013)
#' and Revell (2012).
#'
#' @references
#' Abouheif E. (1999) A method for testing the assumption of phylogenetic independence in comparative data. Evolutionary Ecology Research 1, 895-909.
#' Blomberg S.P., Garland Jr T. & Ives A.R. (2003) Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57, 717-745.
#' Gittleman J.L. & Kot M. (1990) Adaptation: Statistics and a null model for estimating phylogenetic effects. Systematic Biology 39, 227-241.
#' Pagel M. (1999) Inferring the historical patterns of biological evolution. Nature 401, 877-884.
#' Revell L.J. (2012) phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3, 217-223.
#' Pavoine S. & Ricotta C. (2013) Testing for Phylogenetic Signal in Biological Traits: The Ubiquity of Cross-Product Statistics. Evolution 67, 828-840.
#'
#' @seealso \code{\link{phyloSimSignal}} .
#' 
#' @examples
#' require(ape)
#' require(phylobase)
#' data(navic)
#' tipData(navic)$rand <- rnorm(17)
#' tipData(navic)$BM <- rTraitCont(as(navic, "phylo"))
#' phyloSignal(navic)
#' 
#' @rdname phylosignalStats
#' @export
phyloSignal <- function(p4d, methods = c("all", "I", "Cmean", "Lambda", "K", "K.star"), reps = 999, W = NULL){
  methods <- match.arg(methods, several.ok = TRUE)
  p4 <- extractTree(p4d)
  phy <- as(p4, "phylo")
  X <- tdata(p4d, type = "tip")
  X <- as.matrix(X)
  res <- list()
  
  if("all" %in% methods){
    methods <- c("I", "Cmean", "Lambda", "K", "K.star")
  }
  
  if("Cmean" %in% methods){
    WA <- proxTips(p4d, method = "Abouheif", useC = TRUE)
    tmp <- apply(X, 2, moranTest, Wr = WA, reps = reps)
    res$stat$Cmean <- unlist(sapply(tmp, "[", 1))
    res$pvalue$Cmean <- unlist(sapply(tmp, "[", 2))
  }
  
  if("I" %in% methods){
    if(is.null(W)){
      W <- proxTips(p4d, method = "patristic", useC = TRUE)
    } else {
      if(is.matrix(W)){
        W <- W[rownames(X), rownames(X)]
      }
    }
    tmp <- apply(X, 2, moranTest, Wr = W, reps = reps)
    res$stat$I <- unlist(sapply(tmp, "[", 1))
    res$pvalue$I <- unlist(sapply(tmp, "[", 2))
  }
  
  if(any(c("K", "K.star", "Lambda") %in% methods)){  
    VCV <- vcv.phylo(phy, model = "Brownian")
    
    if("K" %in% methods){
      tmp <- apply(X, 2, kTest, vcv = VCV, reps = reps)
      res$stat$K <- unlist(sapply(tmp, "[", 1))
      res$pvalue$K <- unlist(sapply(tmp, "[", 2))
    }
    
    if("K.star" %in% methods){
      tmp <- apply(X, 2, kStarTest, vcv = VCV, reps = reps)
      res$stat$K.star <- unlist(sapply(tmp, "[", 1))
      res$pvalue$K.star <- unlist(sapply(tmp, "[", 2))
    }
    
    if("Lambda" %in% methods){
      tmp <- apply(X, 2, lambdaTest, vcv = VCV)
      res$stat$Lambda <- unlist(sapply(tmp, "[", 1))
      res$pvalue$Lambda <- unlist(sapply(tmp, "[", 2))
      res$pvalue$Lambda <- ifelse(res$pvalue$Lambda < 1/(reps + 1), 1/(reps + 1), res$pvalue$Lambda)
    }
  }
  
  res$stat <- as.data.frame(res$stat, row.names = colnames(X))
  res$pvalue <- as.data.frame(res$pvalue, row.names = colnames(X))
  return(res)
}


#' Test Pagel's Lambda
#' Optimize Pagel's Lambda and do a likelihood ratio test.
#' 
#' @param x a vector of numeric data.
#' @param vcv the phylogenetic variance-covariance matrix.
#' 
#' @details The optimization process is currently performed in R.
#' Could be interesting to do this in Cpp and merge with code{pagelLogLik}.
lambdaTest <- function(x, vcv){
  lambda.max <- max(vcv)/max(vcv[lower.tri(vcv)])
  opt <- suppressWarnings(
    tryCatch(
      optimize(pagelLogLik, interval=c(0, lambda.max), xr = x, vcvr = vcv, maximum = TRUE),
             warning = function(w) return(list(maximum = NA, objective = NA))
      )
    )
  Lambda <- opt$maximum
  if(is.finite(opt$objective)){
    logL0 <- pagelLogLik(0, x, vcv)
    if(is.finite(logL0)){
      pvalue <- as.numeric(pchisq(2 * (opt$objective - logL0), df = 1, lower.tail = FALSE))
    } else {
      pvalue <- NA
    }
  } else {
    pvalue <- NA
  }
  return(list(Lambda = Lambda, pvalue = pvalue))
}


