HWAlltests <- function (x, verbose = TRUE, include.permutation.test = FALSE, x.linked = FALSE) 
{
  # Does all tests for HWE and summarizes the results.    
  out.nc <- HWChisq(x, cc = 0, x.linked = x.linked, verbose = FALSE)
  out.cc <- HWChisq(x, cc = 0.5, x.linked = x.linked, verbose = FALSE)
  out.lr <- HWLratio(x, x.linked = x.linked, verbose = FALSE)
  out.ex.sel <- HWExact(x, pvaluetype = "selome", x.linked = x.linked, verbose = FALSE)
  if(!x.linked) out.ex.dos <- HWExact(x, pvaluetype = "dost", verbose = FALSE) else {
    out.ex.dos <- list(pval=NA)
  }
  out.ex.mid <- HWExact(x, pvaluetype = "midp", x.linked = x.linked, verbose = FALSE)
  if (include.permutation.test) {
    set.seed(123)
    out.perm <- HWPerm(x, verbose = FALSE, x.linked = x.linked)
  }
  else {
    out.perm <- list(stat = NA, pval = NA)
  }
  chisqpval <- out.nc$pval
  chisqccpval <- out.cc$pval
  teststring <- c("Chi-square test:", "Chi-square test with continuity correction:", 
                  "Likelihood-ratio test:", "Exact test with selome p-value:", 
                  "Exact test with dost p-value:", "Exact test with mid p-value:", 
                  "Permutation test:")
  pvals <- c(out.nc$pval, out.cc$pval, out.lr$pval, out.ex.sel$pval, 
             out.ex.dos$pval, out.ex.mid$pval, out.perm$pval)
  stats <- c(out.nc$chisq, out.cc$chisq, out.lr$G2, NA, NA, 
             NA, out.perm$stat)
  Res <- data.frame(stats, pvals)
  colnames(Res) <- c("Statistic", "p-value")
  rownames(Res) <- teststring
  if (verbose) 
    print(Res)
  out <- Res
}
