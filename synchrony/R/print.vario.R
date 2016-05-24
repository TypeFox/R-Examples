print.vario <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {

  switch(x$metric,
         semivar={name='variogram (semivariance)'},
         cov={name='covariogram (covariance)'},
         pearson={name='correlogram (Pearson)'},
         spearman={name='correlogram (Spearman)'},
         kendall={name='correlogram (Kendall)'},
         moran={name='correlogram (Moran\'s I)'},
         geary={name='correlogram (Geary\'s C)'})
  
  if(x$is.multivar)
    mult.or.univ='Multivariate'
  else
    mult.or.univ='Univariate'
  
  if (x$is.centered)
    cent="centered"
  else
    cent="uncentered"
  
  cat(paste(mult.or.univ, cent, name))
  cat("\nMean bin distance:\n")
  print(x$mean.bin.dist, digits=digits)
  cat("\nNumber of points per bin:\n")
  print(x$npoints, digits=digits)
  cat("\nVariogram values:\n")
  print(x$vario, digits=digits)
  
  if (!is.null(x$rands)) {
    switch(x$mult.test.corr,
           none={pvals='Non-adjusted p-values'},
           bonferroni={pvals='Bonferroni adjusted p-values'},
           holm={pvals='Holm-Bonferroni adjusted p-values'},
           hochberg={pvals='Hochberg-adjusted p-values'}
    )
    
    if (x$alternative=='one.tailed')
      tail="one-tailed test"
    else
      tail="two-tailed test"
    cat(paste0("\n", pvals, " based on ",  NROW(x$rands)-1, " randomizations (", tail, "):\n"))
    print(x$pvals, digits=digits)
  }
   
}
