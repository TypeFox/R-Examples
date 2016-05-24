#  print.PMCMR.R
#
#  Copyright (C) 2015, 2016 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

summary.PMCMR <-
function(object, ...)
{
    OK <- inherits(object, "PMCMR")
    if (!OK)
        stop ("Not an object of class PMCMR")
    if (!is.matrix(object$statistic))
        stop ("Matrix object$statistic not found.")
    pval <- as.numeric(object$p.value)
    stat <- as.numeric(object$statistic)
    grp1 <- as.numeric(c(col(object$p.value)))
    cnam <- colnames(object$p.value)
    grp2 <- as.numeric(c(row(object$p.value)))
    rnam <- rownames(object$p.value)
    H0 <- paste(cnam[grp1], " = ", rnam[grp2])
    OK <- !is.na(pval)
    cat("\n\tPairwise comparisons using", object$method, "\n\n")
    cat("data: ", object$data.name, "\n\n")
    cat("\nP value adjustment method:", object$p.adjust.method, "\n")
    xdf <- data.frame(H0 = H0[OK], statistic = stat[OK],
                      p.value = format.pval(pval[OK], 2))
#    out <- list(x = xdf,
#                p.adjust.method = object$p.adjust.method,
#                method = object$method,
#                data.name = object$data.name)
#    class(out) <- "summary.PMCMR"
#    out
    print(xdf)
    invisible(object)
}

#print.summary.PMCMR <-
#function(x, ...) {
#    cat("\n\tPairwise comparisons using", x$method, "\n\n")
#    cat("data: ", x$data.name, "\n\n")
#    cat("\nP value adjustment method:", x$p.adjust.method, "\n")
#    print(x$xdf)
#    invisible(x)
#}

## This was taken from package stat
## file pairwise.R
## (C) 2014 R Core Team, GPL >= 2
print.PMCMR <-
function(x, ...)
{
    cat("\n\tPairwise comparisons using", x$method, "\n\n")
    cat("data: ", x$data.name, "\n\n")
    pp <- format.pval(x$p.value, 2, na.form="-")
    attributes(pp) <- attributes(x$p.value)
    print(pp, quote=FALSE, ...)
    cat("\nP value adjustment method:", x$p.adjust.method, "\n")
    invisible(x)
}

get.pvalues <-
function(object, ...)
{
    OK <- inherits(object, c("PMCMR", "pairwise.htest"))
    if (!OK)
        stop ("Not an object of class PMCMR or pairwise.htest")
    if (!is.matrix(object$p.value))
        stop ("Matrix object$p.value not found.")
    pval <- as.numeric(object$p.value)
    grp1 <- as.numeric(c(col(object$p.value)))
    cnam <- colnames(object$p.value)
    grp2 <- as.numeric(c(row(object$p.value)))
    rnam <- rownames(object$p.value)
    H0 <- paste(cnam[grp1],"-",rnam[grp2], sep="")
    OK <- !is.na(pval)
    out <- pval[OK]
    names(out) <- H0[OK]
    out
}
