
uroot.raw.pvalue <- function(x, type = c("CH", "HEGY"), v, n, ctd, S, Ftpi)
{
  type <- match.arg(type)

  switch(type,

    "CH" = 
    {
      if (v > 12)
        stop("tabulated values are not available for ", sQuote("v"), " = ", v, " degrees of freedom")
      pval <- approx(c(0, .CH.orig.cv[v,], Inf), 
        c(1, 0.20, 0.10, 0.075, 0.05, 0.025, 0.01, 0), x)$y
      return(pval)
    },

    "HEGY" = 
    {
      if (S == 4) {
        tab <- .HEGY.orig.cv[ctd,,,Ftpi]
        nsizes <- c(48, 100, 136, 200)
      } else 
      if (S == 12) {
        tab <- .BM.orig.cv[ctd,,,Ftpi]
        nsizes <- c(240, 480, 10000)
      } else {
        return(NA)
        warning("tabulated values are available only for quarterly and monthly series. ") 
      }

      #option rule=2 chooses the values tabulated for n=48 and n=200 if "n" is beyond these limits
      cvals <- rep(NA, 4)
      for (i in seq_along(cvals))
        cvals[i] <- approx(nsizes, c(tab[,i]), n, rule = 2)$y
      if (Ftpi == "pair") {
        pval <- approx(c(10000, cvals, 0), c(0, 0.01, 0.025, 0.05, 0.1, 1), x)$y
      } else 
        pval <- approx(c(-10000, cvals, 3), c(0, 0.01, 0.025, 0.05, 0.1, 1), x)$y
      return(pval)
    }
  ) # end switch
}
