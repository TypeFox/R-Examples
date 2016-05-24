#' S3 class solarAssoc.
#'
#' @name solarAssocClass
#' @rdname solarAssocClass
#'
#' @param x 
#'    An object of class \code{solarAssoc}.
#' @param object
#'    An object of class \code{solarAssoc}.
#' @param y
#'    Character argument for \code{plot} method, 
#'    indicating the type of plot.
#'    The default value is \code{"manh"}.
#' @param alpha
#'    Numeric argument between 0 and 1 for \code{summary} method, 
#'    indicating the significance level after Bonferroni multiple-test correction.
#'    The default value is 0.05.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass solarAssoc

#--------------------
# Print method
#--------------------

#' @rdname solarAssocClass
#' @export
print.solarAssoc <- function(x, ...)
{
  cat("\nCall: ")
  print(x$assoc$call)
  
  cat("\n Input SNP data:\n")
  switch(x$assoc$assoc.informat,
    "snpdata" = cat("  *  ", x$assoc$num.snps, " SNP genotypes passed by `snpdata` argument\n", sep = ""),
    "snpcovdata" = cat("  *  ", x$assoc$num.snps, " SNP covariates passed by `snpcovdata` argument\n", sep = ""),
    "genocov.file" = cat("  *  SNP covariates passed in ", length(x$assoc$genocov.files),
      " file(s) by `genocov.files` argument\n", sep = ""),
    "genocov.files" = cat("  *  SNP covariates passed in ", length(x$assoc$genocov.files),
      " file(s) by `genocov.files` argument and ", length(x$assoc$snplists.files), 
      " files(s) by `snplists.files` argument\n", sep = ""),
    stop("switch error")
  )
  
  cat("\n Output results of association:\n")
  if(x$assoc$assoc.outformat == "df") {
    cat("\n  *  Table of association results (first 5 out of ", nrow(x$snpf), " rows):\n", sep = "")
    print(head(x$snpf, 5))
  }
  #cat("  *  assoc.outformat: ", x$assoc$assoc.outformat, "\n", sep = "")
  
  batches.str <- ""
  nb <- modelParNumBatches(x)
  if(nb > 1) {
    batches.str <- paste0(", ", nb, " batches")
  }
  cat("\n CPU time on ", modelParCores(x), " core(s)",
    batches.str,
    ": ", 
    modelParCPUtime(x, "POSIX"), "\n", sep = "")
}

#' @rdname solarAssocClass
#' @export
plot.solarAssoc <- function(x, y = "manh", ...)
{
  switch(y,
    manh = plotManh(x, ...),
    qq = plotQQ(x, ...),
    stop("switch error"))
}

plot.solarAssoc.old <- function(x, 
  alpha = 0.05, corr = "BF", pval = "pval", ...)
{
  #ret <- require(ggplot2)
  #if(!ret) {
  #  stop("`ggplot2` package is required for plotting")
  #}
  ret <- requireNamespace("scales", quietly = TRUE)
  if(!ret) {
    stop("`scales` package is required for plotting")
  }
    
  df <- x$snpf
  N <- nrow(df)
  
  ord <- order(df$pSNP)
  df <- df[ord, ]
  df$ord <- 1:nrow(df)
    
  ### multiple-test correction  
  pSNP <- NULL # R CMD check: no visible binding
  df <- within(df, {
    qSNP <- switch(corr,
      "BF" = pSNP * nrow(df),
      stop("switch error (1) in `plot.solarAssoc`"))
  })

  df <- mutate(df,
    signif = qSNP < alpha)
  
  num.signif <- sum(df$signif)
  
  ### subset
  num.nonsignif <- ifelse(num.signif, 3, 5)
  ord.nonsignif <- seq(1, num.nonsignif) + num.signif
  
  ### `pf` data frame for plotting
  pf <- rbind(subset(df, signif),
    subset(df, ord %in% ord.nonsignif))
  
  ord <- order(pf$ord)
  pf <- pf[ord, ]
  
  ### plotting settings
  # custom axis transformations in ggplot2
  #  - @ http://www.numbertheory.nl/2012/08/14/custom-axis-transformations-in-ggplot2/
  #  - @ http://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
  reverselog_trans <- function(base = exp(1)) 
  {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
      scales::log_breaks(base = base), 
      domain = c(1e-100, Inf))
  }
  
  xbreaks <- pf$ord
  xlables <- pf$SNP

  p <- switch(pval,
    "pval" = ggplot(pf, aes(ord, pSNP)) + geom_point() + 
      geom_segment(aes(x = ord, xend = ord, y = 1, yend = pSNP)) + 
      geom_hline(yintercept = alpha/N),
    "qval" = ggplot(pf, aes(ord, qSNP)) + geom_point() + 
      geom_segment(aes(x = ord, xend = ord, y = 1, yend = qSNP)) + 
      geom_hline(yintercept = alpha/N),
      stop("switch error (2) in `plot.solarAssoc`"))

  ylab <- switch(pval,
    "pval" = paste("P-value (alpha = ", format(alpha), ", ", "alpha/N = ", format(alpha/N), ")", sep = ""),
    "qval" = paste("Q-value (corrected p-value) (alpha = ", format(alpha), ")", sep = ""),
     stop("switch error (3) in `plot.solarAssoc`"))
  
  title <- paste("Association model: ", num.signif, " significant, N = ", N, sep = "")
    
  p <- p + scale_y_continuous(trans = reverselog_trans(10)) + 
    scale_x_continuous(breaks = xbreaks, labels = xlables) + 
    labs(y = ylab, x = "SNP", title = title) +
    coord_flip()
  
  return(p)
}

#' @rdname solarAssocClass
#' @export
summary.solarAssoc <- function(object, alpha = 0.05, ...)
{
  cat("\nCall: ")
  print(object$assoc$call)
  
  ### var
  num.snps <- nrow(object$snpf)
  
  cat("\nAssociation model\n")
  cat(" * Number of SNPs:", num.snps, "\n")
  cat(" * Input format:", object$assoc$assoc.informat, "\n")
  
  # signif. SNPs
  pSNP.thr <- alpha / num.snps
  pSNP <- NULL # R CMD check: no visible binding
  snpf <- subset(object$snpf, pSNP < pSNP.thr)
  num.snps.signif <- length(which(object$snpf$pSNP<pSNP.thr))
  cat(" * Number of significal SNPs: ", num.snps.signif, 
    " (Bonferroni correction with alpha ", alpha, ")\n", sep = "")
  if(num.snps.signif > 0) {
    ord <- with(snpf, order(pSNP))
    snpf <- snpf[ord, ]
    print(snpf, 5)
  }
}

#--------------------
# Other methods
#--------------------

#' Annotate SNPs in association study
#'
#' The function calls \code{\link{annotateSNPs}} function,
#' which does the job.
#'
#' @seealso \code{\link{annotateSNPs}}
#'
#' @param x 
#'    An object of class \code{solarAssoc} or a character vector of SNPs.
#' @param mode
#'    A character with the mode of SNPs selection.
#'    Possible values are \code{"significant"}, \code{"top"} and \code{"all"}.
#'    The default value is \code{"significant"}.
#' @param alpha
#'    A numeric value from 0 to 1, the significance level after Bonferroni multiple-test correction.
#'    Corresponds to \code{mode} equal to \code{"significant"}.
#' @param num.top
#'    An integer value, the number of top SNPs to be annotated.
#'    Corresponds to \code{mode} equal to \code{"top"}.
#'    The default value is 10.
#' @param ...
#'    Additional arguments passed to \code{annotateSNPs}.
#' @return
#'    A data table with annotation results.
#'
#' @export
annotate <- function(x, mode = c("significant", "top", "all"), 
  alpha = 0.05,
  num.top = 10, ...)
{
  annot <- annotateSNPs(x, mode = mode, alpha = alpha, num.top = num.top, ...)
#      Query Chromosome    Marker Class Gene Alleles Major Minor   MAF        BP
#1 rs2731672          5 rs2731672   snp  F12     A/G     A     G 0.484 177415472

  ### Case 1: annotate the association results in `x`
  ### merge
  #annot <- data.table(annot)
  #setnames(annot, "Query", "SNP")

  #setkey(A$snpf, SNP)
  #setkey(annot, SNP)

  #A$snpf <- merge(A$snpf, annot, all.x = TRUE)

  ### Case 2: annotate the selected SNPs
  annot <- data.table(annot)
  setnames(annot, "Query", "SNP")
  
  SNP <- NULL # fix `no visible binding`
  setkey(x$snpf, SNP)
  setkey(annot, SNP)

  annot <- merge(x$snpf, annot, all.x = FALSE, all.y = TRUE) # names of `x$snpf` go first

  return(annot)  
}
