## plot.genotype.R
###------------------------------------------------------------------------
## What: Plot genotype object
## $Id: plot.genotype.R 1274 2007-07-18 12:09:37Z ggorjan $
## Time-stamp: <2007-07-18 16:06:07 ggorjan>
###------------------------------------------------------------------------

plot.genotype <- function(x,
                          type=c("genotype", "allele"),
                          what=c("percentage","number"),
                          ...)
{
  what <- match.arg(what)
  type <- match.arg(type)

  ## get details
  tmp <- summary(x)

  ## Percentages or numbers
  whati <- ifelse(what == "percentage", 2, 1)

  ## Plot
  if (type == "allele") {
    barplot(tmp$allele.freq[, whati], ...)
  } else { # genotype
    barplot(tmp$genotype.freq[, whati], ...)
  }
}

###------------------------------------------------------------------------
## plot.genotype.R ends here
