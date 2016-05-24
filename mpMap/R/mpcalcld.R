#' Calculate linkage disequilibrium between all pairs of markers
#'
#' Calculate linkage disequilibrium for multiple multi-allelic measures, including Lewontin's D, W, and correlation. 
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @return Original object with additional component ld consisting of separate matrices for each different measure.
#' \item{W}{Related to Chi-square test for independence (measure of non-centrality)}
#' \item{LewontinD}{Lewontin's D' statistic}
#' \item{delta2}{Delta-squared statistic}
#' \item{r2}{Correlation between markers} 
#' @seealso \code{\link[mpMap]{mpestrf}}
#' @references Zhao H, Nettleton D, Soller M and JCM Dekkers (2005) Evaluation of linkage disequilibrium measures between multi-allelic markers as predictors of linkage disequilibrium between markers and QTL. Genet Res Camb 86:77-87.

mpcalcld <-
function(object)
{ 
  if (missing(object))	
	stop("Missing a required argument for this function")

  if (!inherits(object, "mpcross")) stop("Object must be of class mpcross")

  if (is.null(object$rf)) stop("Must compute recombination fractions first")

  require(gdata)
 
  output <- object 
  final.g <- cbind(object$id, object$finals)
  rmat <- object$rf$theta
  rmat[is.na(rmat)] <- 0.5

  n.loci <- ncol(object$founders)

  # pairs of loci
  if (n.loci>1) {
    pairsij <- t(cbind(combn(1:n.loci, 2), matrix(data=rep(1:n.loci,2), nrow=2, ncol=n.loci, byrow=TRUE)))
  pairsij <- pairsij[order(pairsij[,1], pairsij[,2]), c(2:1)]
  } else pairsij <- matrix(data=rep(1:n.loci,each=2), nrow=n.loci, ncol=2, byrow=TRUE)

  ld <- CR_calcLD(final.g, object$founders, object$pedigree, pairsij, rmat)

  output$ld <- list()
  output$ld$W <- matrix(nrow=n.loci, ncol=n.loci) 
  output$ld$LewontinD <- matrix(nrow=n.loci, ncol=n.loci)
  output$ld$delta2 <- matrix(nrow=n.loci, ncol=n.loci)
#  output$ld$r2 <- cor(object$finals, use="pairwise.complete.obs")^2
  output$ld$r2 <- matrix(nrow=n.loci, ncol=n.loci)

  lowerTriangle(output$ld$W, diag=TRUE) <- ld[,1]
  lowerTriangle(output$ld$LewontinD, diag=TRUE) <- ld[,2]
  lowerTriangle(output$ld$delta2, diag=TRUE) <- ld[,3]
  lowerTriangle(output$ld$r2, diag=TRUE) <- ld[,4]

  upperTriangle(output$ld$W) <- upperTriangle(t(output$ld$W))
  upperTriangle(output$ld$LewontinD) <- upperTriangle(t(output$ld$LewontinD))
  upperTriangle(output$ld$delta2) <- upperTriangle(t(output$ld$delta2))
  upperTriangle(output$ld$r2) <- upperTriangle(t(output$ld$r2))

  output$ld$W[is.na(object$rf$theta)] <- NA
  output$ld$LewontinD[is.na(object$rf$theta)] <- NA
  output$ld$delta2[is.na(object$rf$theta)] <- NA
  output$ld$r2[is.na(object$rf$theta)] <- NA

  rownames(output$ld$W) <- rownames(output$ld$LewontinD) <- rownames(output$ld$delta2) <- rownames(output$ld$r2) <- colnames(output$ld$W) <- colnames(output$ld$LewontinD) <- colnames(output$ld$delta2) <- colnames(output$ld$r2) <- colnames(object$finals)

  return(output)
}

