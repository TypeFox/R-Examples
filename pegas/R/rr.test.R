## rr.test.R (2011-07-11)

##   Tajima Relative Rate Test of Molecular Clock

## Copyright 2011 Alastair Potts

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

rr.test <- function(x, y, out)
{
    l <- seg.sites(rbind(x, out))
    m <- seg.sites(rbind(y, out))
    m1 <- length(l[! l %in% m])
    m2 <- length(m[! m %in% l])
    CHI <- (m1 - m2)^2 / (m1 + m2)
    PVAL <- 1 - pchisq(CHI, df = 1)
    list(Chi = CHI, Pval = PVAL)
}
