# $Id: HWE.chisq.R 1352 2012-08-14 14:21:35Z warnes $


###
### Hardy-Weinberg Equilibrium Significance Test
###


HWE.chisq <- function(x, ...)
  UseMethod("HWE.chisq")

HWE.chisq.genotype <- function (x, simulate.p.value = TRUE, B = 10000, ...)
{
    observed.no <- table(factor(allele(x, 1), levels = allele.names(x)), 
        factor(allele(x, 2), levels = allele.names(x)))
    tab <- observed.no
    tab <- 0.5 * (tab + t(tab))
    k <- ncol(tab)
    if(simulate.p.value)
      {
        test <- chisq.test(tab, simulate.p.value = simulate.p.value, 
                           B = B, ...)
      }
    else
      {
        test <- chisq.test(tab, ...)
        test$parameter <- k*(k-1)/2
        test$p.value <- pchisq(test$statistic, test$parameter, lower.tail = FALSE)
        names(test$statistic) <- "X-squared"
        names(test$parameter) <- "df"
      }
    return(test)
}

