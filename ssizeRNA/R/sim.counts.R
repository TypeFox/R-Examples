#' RNA-seq Count Data Simulation from Negative-Binomial Distribution
#' 
#' This function simulates count data from Negative-Binomial distribution
#' for two-sample RNA-seq experiments with given mean, dispersion 
#' and log fold change. 
#' A count data matrix is generated.
#' 
#' @import MASS
#' @param nGenes total number of genes, the default value is \code{10000}.
#' @param pi0 proportion of non-differentially expressed genes, 
#'            the default value is \code{0.8}.
#' @param m sample size per treatment group.
#' @param mu a vector (or scalar) of mean counts in control group 
#'           from which to simulate.
#' @param disp a vector (or scalar) of dispersion parameter 
#'             from which to simulate.
#' @param logfc a vector (or scalar, or a function that takes an integer n 
#'                        and generates a vector of length n)
#'              of log fold change for differentially expressed (DE) genes.  
#' @param up proportion of up-regulated genes among all DE genes, 
#'           the default value is \code{0.5}.
#' @param replace sample with or without replacement from given parameters. 
#'                See Details for more information.
#' 
#' @details If the total number of genes \code{nGenes} is larger 
#'          than length of \code{mu} or \code{disp}, 
#'          \code{replace} always equals \code{TRUE}.
#' 
#' @return \item{counts}{RNA-seq count data matrix.}
#' @return \item{group}{treatment group vector.}
#' @return \item{lambda0}{mean counts in control group for each gene.}
#' @return \item{phi0}{dispersion parameter for each gene.}
#' @return \item{de}{differentially expressed genes indicator: 
#'                   \code{0} for non-differentially expressed genes, 
#'                   \code{1} for up-regulated genes, 
#'                   \code{-1} for down-regulated genes.}
#' @return \item{delta}{log fold change for each gene between 
#'                      treatment group and control group.}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, 
#'         Peng Liu \email{pliu@@iastate.edu}
#' 
#' @examples
#' m <- 3                    ## sample size per treatment group
#' mu <- 10                  ## mean counts in control group for all genes 
#' disp <- 0.1               ## dispersion for all genes
#' logfc <- log(2)           ## log fold change for DE genes
#' 
#' sim <- sim.counts(m = m, mu = mu, disp = disp, logfc = logfc)
#' sim$counts                ## count data matrix
#' 
#' 
#' ## varying fold change
#' logfc1 <- function(x){rnorm(x, log(2), 0.5*log(2))}
#' sim1 <- sim.counts(m = m, mu = mu, disp = disp, logfc = logfc1)
#' 
#' @export
#' 
sim.counts <- function(nGenes = 10000, pi0 = 0.8, m, mu, disp, logfc, 
                       up = 0.5, replace = TRUE){
  arg <- list(nGenes = nGenes,
              pi0 = pi0,
              group = rep(c(1, 2), each = m))
 
  ## expected false positives
  FP <- round(nGenes * pi0)
  TP <- nGenes - FP 
  
  ## types of true positives
  TP_up <- round(TP * up)
  TP_down <- TP - TP_up 

  de <- c(rep(0, FP), rep(1, TP_up), rep(-1, TP_down))
  de <- de[sample.int(length(de))] ## resample
  
  # h = vector indicating which pseudo-genes to re-simulate
  h <- rep(TRUE, nGenes) 
  counts <- matrix(0, nrow = nGenes, ncol = 2 * m)
  
  ## log fold change, approximately half positive, half negative
  delta <- rep(0, nGenes)
  if (is.function(logfc)){
    lfc <- logfc(TP)
  }else{
    lfc <- logfc
  }
  delta[de != 0] <- lfc * de[de != 0]
  
  selected_genes <- true_means <- true_disps <- rep(0, nGenes)
  left_genes <- 1:length(mu)
  lambda <- phi <- matrix(0, nrow = nGenes, ncol = 2 * m)
  
  while(any(h)){
    temp <- sample.int(length(left_genes), sum(h), replace)
    temp <- temp[order(temp)]
    selected_genes[h] <- left_genes[temp]
    if (replace == FALSE){
      left_genes <- left_genes[-temp]
    }
    
    true_means[h] <- mu[selected_genes[h]]
    true_disps[h] <- disp[selected_genes[h]]
    
    lambda[h,] <- matrix(true_means[h], ncol = 1) %*% 
                  matrix(rep(1, 2 * m), nrow = 1) * 
                  cbind(matrix(rep(1, sum(h) * m), ncol = m), 
                        matrix(rep(exp(delta[h]), m), ncol = m))
    ## mean of counts
    
    phi[h,] <- matrix(rep(true_disps[h], 2 * m), ncol = 2 * m)
    ## dispersion of counts
    
    counts[h,] <- rnegbin(sum(h) * 2 * m, lambda[h,], 1 / phi[h,])
    h <- (rowSums(cpm(counts) > 2) < 3)
    # print(sum(h))
  }
  
  if(any(rowSums(cpm(counts) > 2) < 3 ))
    print("Error: Failed to simulate data: some genes are not expressed.")
  
  list(counts = counts, 
       group = arg$group, 
       lambda0 = lambda[, 1],   # mean counts in control group
       phi0 = phi[, 1],   # dispersion
       de = de,   # DE indicator
       delta = delta  # log fold change
  )
}
