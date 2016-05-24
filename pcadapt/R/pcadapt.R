#' Principal Component Analysis for outlier detection
#'
#' \code{pcadapt} performs principal component analysis and computes p-values to test for outliers. The test for
#' outliers is based on the correlations between genetic variation and the first \code{K} principal components.
#' \code{pcadapt} also handles Pool-seq data for which the statistical analysis is
#' performed on the genetic markers frequencies. Returns an object of class \code{pcadapt}.
#'
#' @details First, a principal component analysis is performed on the scaled and centered genotype data. To account for missing
#' data, the correlation matrix between individuals is computed using only the markers available for each
#' pair of individuals. Depending on the specified \code{method}, different test statistics can be used.
#'
#' \code{mahalanobis} (default): the Mahalanobis distance is computed for each genetic marker using a robust
#' estimate of both mean and covariance matrix between the \code{K} vectors of z-scores.
#'
#' \code{communality}: the communality statistic measures the proportion of variance explained by the first \code{K} PCs.
#'
#' \code{euclidean}: the Euclidean distance between the \code{K} z-scores of each genetic marker and the mean of the \code{K} vectors of z-scores is computed.
#'
#' \code{componentwise}: returns a matrix of z-scores.
#'
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic inflation factor (\code{gif}) when \code{method="mahalanobis","euclidean"}.
#' When \code{method="communality"}, the test statistic is first multiplied by \code{K} and divided by the percentage of variance explained by the first \code{K} PCs
#' before accounting for genomic inflation factor. When using \code{method="mahalanobis","communality","euclidean"}, the scaled statistics (\code{chi2_stat}) should follow
#' a chi-squared distribution with \code{K} degrees of freedom. When using \code{method="componentwise"}, the z-scores should follow a chi-squared distribution with \code{1}
#' degree of freedom. For Pool-seq data, \code{pcadapt} provides p-values based on the Mahalanobis distance for each SNP.
#'
#' @param input a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Four statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, \code{"euclidean"} and \code{"componentwise"}.
#' @param data.type a character string specifying the type of data being read, either a \code{genotype} matrix (\code{data.type="genotype"}),
#' or a matrix of allele frequencies (\code{data.type="pool"}).
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param transpose a logical value indicating whether the genotype matrix has to be tranposed or not. A genotype
#' matrix should be \code{p x n} where \code{p} is the number of genetic markers and \code{n} is the number of individuals.
#' If the data contains missing values, please encode missing values as \code{9} or use the function \code{read4pcadapt}
#' to format the data.
#' @param output.filename a character string specifying the names of the files created by \code{pcadapt}.
#' @param clean.files a logical value indicating whether the auxiliary files should be deleted or not.
#'
#' @return  The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @importFrom robust covRob
#' @importFrom MASS cov.rob
#' @importFrom stats median pchisq na.omit qchisq
#' @importFrom utils read.table write.table
#'
#' @useDynLib pcadapt wrapper_pcadapt wrapper_converter
#'
#' @export
#'
pcadapt <- function(input,
                        K=2,
                        method="mahalanobis",
                        data.type="genotype",
                        min.maf=0.05,
                        ploidy=2,
                        output.filename="pcadapt_output",
                        clean.files = TRUE,
                        transpose = FALSE){

  #############################################
  ########## test arguments and init ##########
  #############################################
  
  if (data.type == "genotype"){
    
    input.filename <- 1
    
    # input file
    if ((class(input) == "character") && (!file.exists(input))){
      stop(paste0("File ",input," does not exist."))
    } else if ((class(input) == "character") && (file.exists(input))){
      if (transpose == TRUE){
        .C("wrapper_converter",as.character(input),as.integer(2),PACKAGE = "pcadapt")
        input.filename <- paste0(input,".pcadapt")
      } else {
        input.filename <- input
      }
    }
    
    if (class(input) %in% c("array","matrix","data.frame")){
      if (!(ncol(input)>0) || !(nrow(input)>0)){
        stop("Invalid input genotype matrix.")
      }
      if (transpose == TRUE){
        write.table(t(input),"tmp.pcadapt",col.names = FALSE,row.names = FALSE)
      } else {
        write.table(input,"tmp.pcadapt",col.names = FALSE,row.names = FALSE)
      }
      input.filename <- "tmp.pcadapt"
    }
    
    if (class(input.filename) != "character" || (!file.exists(input.filename))){
      stop("Invalid argument. Make sure the file exists or that the data is in the workspace.")
    }
    
    if (class(K) != "numeric" || K <= 0){
      stop("K has to be a positive integer.")
    }
    
    if (!(method %in% c("mahalanobis","communality","euclidean","componentwise"))){
      warning("Unknown method. 'mahalanobis' will be used hence.")
      method <- "mahalanobis"
    }
    
    if (class(min.maf) != "numeric" || min.maf < 0 || min.maf > 0.45){
      warning("min.maf has to be a real number between 0 and 0.45. Default value will be used hence.")
      min.maf <- 0.05
    }
    
    if (!(ploidy %in% c(1,2))){
      stop("pcadapt only supports haploid and diploid data.")
    }
    
    if (class(output.filename) != "character"){
      warning("output.filename has to be a character string. Default value will be used hence.")
      output.filename <- "pcadapt_output"
    }
    
    if (ploidy == 2){
      .C("wrapper_pcadapt",
         as.character(input.filename),
         as.integer(K),
         as.double(min.maf),
         as.integer(0),
         as.character(output.filename),
         PACKAGE = "pcadapt"
      );
    } else if (ploidy == 1){
      .C("wrapper_pcadapt",
         as.character(input.filename),
         as.integer(K),
         as.double(min.maf),
         as.integer(1),
         as.character(output.filename),
         PACKAGE = "pcadapt"
      );
    }
    
    res <- create.pcadapt(output.filename,K,method,data.type,min.maf)
    
    if (clean.files == TRUE){
      # Clean files
      # if (file.exists(input.filename)){
      #   file.remove(input.filename)
      # }
      file.remove(paste0(output.filename,".loadings"))
      file.remove(paste0(output.filename,".scores"))
      file.remove(paste0(output.filename,".zscores"))
      file.remove(paste0(output.filename,".maf"))
      file.remove(paste0(output.filename,".sigma"))
    }
  } else if (data.type == "pool"){
    if (class(input) %in% c("array","matrix","data.frame")){
      data <- input  
    } else if (class(input) == "character"){
      data <- read.table(input)  
    } else {
      stop("Invalid input argument.")
    }
    if (transpose == TRUE){
      data <- t(input)
    } 
    nPOP <- nrow(data)
    nSNP <- ncol(data)
    if (missing(K)){
      K <- nPOP-1
    }
    res <- corpca(data,K)
    freq <- apply(data,2,FUN=function(x){mean(x,na.rm=TRUE)})
    res$maf <- as.vector(pmin(freq,1-freq))
    res$loadings[res$maf<min.maf] <- NA 
    res$stat <- array(NA,dim=nSNP)
    finite.list <- which(!is.na(apply(abs(res$loadings),1,sum)))
    if (K>1){
      res$stat[finite.list] <- as.vector(robust::covRob(res$loadings,na.action=na.omit,estim="pairwiseGK")$dist)
    } else {
      onedcov <- as.vector(MASS::cov.rob(res$loadings[finite.list,1]))
      res$stat <- (res$zscores[,1]-onedcov$center)^2/onedcov$cov[1]
    }
    res$gif <- median(res$stat,na.rm=TRUE)/qchisq(0.5,df=K)
    res$chi2.stat <- res$stat/res$gif
    # Compute p-values
    res$pvalues <- compute.pval(res$chi2.stat,K,method="mahalanobis")
    class(res) <- 'pcadapt'
    attr(res,"K") <- K
    attr(res,"method") <- method
    attr(res,"data.type") <- data.type
    attr(res,"min.maf") <- min.maf
  }
  return(res)
}




