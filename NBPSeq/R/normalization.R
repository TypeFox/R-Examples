##' \code{estimate.norm.factors} estiamtes normalization factors to
##' account for apparent reduction or increase in relative frequencies
##' of non-differentially expressing genes as a result of compensating
##' the increased or decreased relative frequencies of truly
##' differentially expressing genes.
##'
##' We take gene expression to be indicated by relative frequency of
##' RNA-Seq reads mapped to a gene, relative to library sizes (column
##' sums of the count matrix). Since the relative frequencies sum to 1
##' in each library (one column of the count matrix), the increased
##' relative frequencies of truly over expressed genes in each column
##' must be accompanied by decreased relative frequencies of other
##' genes, even when those others do not truly differentially
##' express. If not accounted for, this may give a false impression of biological relevance
##' (see, e.g., Robinson and Oshlack (2010), for some
##' examples.)  A simple fix is to compute the relative frequencies
##' relative to effective library sizes---library sizes multiplied by
##' normalization factors.
##'
##' @title Estiamte Normalization Factors
##' @export
##' @param counts a matrix of RNA-Seq read counts with rows
##' corresponding to gene features and columns corresponding to
##' independent biological samples.
##' @param lib.sizes  a vector of observed library sizes, usually and
##' by default estimated by column totals.
##' @param method a character string specifying the method for
##' normalization, currenlty, can be NULL or "AH2010". If method=NULL, the
##' normalization factors will have values of 1 (i.e., no
##' normalization is applied); if method="AH2010" (default), the normalization
##' method proposed by Anders and Huber (2010) will be used.
##' @references {Anders, S. and W. Huber (2010): "Differential
##' expression analysis for sequence count data," Genome Biol., 11,
##' R106.
##'
##' Robinson, M. D. and A. Oshlack (2010): "A scaling normalization method for differential expression analysis of RNA-seq data," Genome Biol., 11, R25.}
##' @return a vector of normalization factors.
##' @examples
##' ## Load Arabidopsis data
##' data(arab)
##'
##' ## Estimate normalization factors using the method of Anders and Huber (2010)
##' norm.factors = estimate.norm.factors(arab);
##' print(norm.factors);
estimate.norm.factors = function(counts, lib.sizes=colSums(counts), method="AH2010") {
  if (is.null(method)) {
    ## No normalization is used
    norm.factors = rep(1, dim(counts)[2]);
  } else {

    if (method=="AH2010") {
      ## The method proposed by Anders and Huber Genome Biology 2010,
      ## 11:R106.
      norm.factors = estimate.norm.factors.AH2010(counts, lib.sizes);
    } else {
      stop("Unknown normalization method.");
    }

  }

  norm.factors
}

##' Estimate normalization factors using the method of Anders and
##' Huber (2010).
##'
##' The normalization factor of a column is computed as the median
##' fold change in relative read frequencies between this column and a
##' reference column. The frequencies in the reference column is the
##' geometric mean of read frequencies from all columns in the matrix
##' \code{counts}.
##'
##' The normalization factors are normalized so that they multiply to 1.
##' 
##' Rows with any 0 counts will be ignored (since their geometric
##' means will be 0) when computing the median fold change.
##'
##' See \code{\link{estimate.norm.factors}} for discussion on the need
##' for normalization.
##' 
##' @title (private) Estimate normalization factors using the method
##' of Anders and Huber (2010).
##' @noRd
##' @param counts a matrix, frequencies of RNA-Seq reads
##' @param lib.sizes a vector, estiamted library sizes
##' @return a vector of normalization factors.
estimate.norm.factors.AH2010 = function(counts, lib.sizes=colSums(counts)) {
  ## The method proposed by Anders and Huber Genome Biology 2010,
  ## 11:R106.

  ## Compute relative frequences, same dimension as counts
  r = t(t(counts)/lib.sizes);

  ## Create a reference column (geometric mean of the read
  ## frequenices in each column)
  ref = exp(rowMeans(log(counts)));
  r.ref = ref / sum(ref);

  ## Estimate median fold change relative to the reference column
  norm.factors = apply(r,  2, function(y) median((y/r.ref)[r.ref>0]));

  ## Normalize the normalizaton factors, so that they multiply to 1.
  norm.factors/exp(mean(log(norm.factors)));
}
