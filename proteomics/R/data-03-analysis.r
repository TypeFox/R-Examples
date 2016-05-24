# Copyright (C) 2012-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

######################################################
# Analysing responses by ANOVA on protein level data #
######################################################

#' Data Analysis -- Testing one feature without Tukey Honest Significant Differences
#'
#' @param model ANOVA model of the corresponding fit
#' @export
testingOneshot <- function (model) {
    pvalue <- anova(model)[["Pr(>F)"]]
    pvalue <- pvalue[-length(pvalue)]
    return(data.frame(.id=c("typus", "treatment"), pvalue=pvalue))
}

testingWithTukey <- function (tuk) {
    mat <- as.data.frame(tuk)
    mat$order <- with(mat, order(order(abs(mat$diff), decreasing=TRUE)))
    mat$seen <- rownames(mat)
    rownames(mat) <- NULL
    mat <- rename(mat, c("p adj"="pvalue"))
    return(mat)
}

#' Data Analysis -- Testing one feature with Tukey Honest Significant Differences
#'
#' @param model ANOVA model of the corresponding fit
#' @param conf.level confidence level
#' @export
testingTukey <- function (model, conf.level) {
    tuk <- TukeyHSD(x=model, conf.level=conf.level)
    mat <- ldply(tuk, testingWithTukey)
    return(mat)
}

#' Data Analysis -- Testing features with Tukey Honest Significant Differences
#'
#' @param dp iTRAQ data in long format
#' @param frm formal for the test
#' @param conf.level confidence level
#' @export
testing <- function(dp, frm, conf.level) {
    tt <- testingTukey(aov(frm, data=dp), conf.level)
    tt$accum.spectra <- nlevels(droplevels(dp$id))
    tt$accum.peptides <- nlevels(droplevels(dp$peptide))
    tt$accum.sequences <- nlevels(droplevels(dp$sequence))
    return(tt)
}

#' Data Analysis -- Testing on peptide level
#'
#' @param dat iTRAQ data in long format
#' @param frm formal for the test
#' @param conf.level confidence level
#' @param ... arguments understood by ddply
#' @export
testForPeptideEffect <- function (dat, frm, conf.level, ...) {
    res <- ddply(dat, .(protein, peptide), testing, ..., frm=frm, conf.level=conf.level)
    res <- rename(res, c(.id="factors"))
    res$factors <- with(res, factor(factors))
    res$seen <- with(res, factor(seen))
    return(res)
}

#' Data Analysis -- Testing on protein level
#'
#' @param dat iTRAQ data in long format
#' @param frm formal for the test
#' @param conf.level confidence level
#' @param ... arguments understood by ddply
#' @export
testForProteinEffect <- function (dat, frm, conf.level, ...) {
    res <- ddply(dat, .(protein), testing, ..., frm=frm, conf.level=conf.level)
    res <- rename(res, c(.id="factors"))
    res$factors <- with(res, factor(factors))
    res$seen <- with(res, factor(seen))
    return(res)
}

##############
## Plotting ##
##############

#' Plotting p-value distributions
#'
#' @param restest result frame of test results
#' @export
pAction <- function (restest) {
    tt <- ggtitle(paste(strwrap("Kernel density estimates of p-values testing
                                the null hypothesis of non existinig
                                differences", width=50), collapse="\n"))
    return(  ggplot(restest, aes(x=pvalue))
           + geom_density(aes(fill=factors), alpha=.3)
           + theme(legend.position="none")
           + facet_wrap(~seen)
           + tt
           )
}

#' Volcano plot
#'
#' @param res result frame of test results
#' @param threshold for biological reasonable effect
#' @param .foldchange wheather results given in ratios or log-ratios
#' @param .plot if true adds a plotting layer
#' @export
pVolcano <- function (res, threshold, .foldchange=TRUE, .plot=TRUE) {
    if (.foldchange) threshold <- log2(threshold)
    p <- (  ggplot(res, aes(diff, -10*log10(pvalue)))
          + geom_vline(xintercept=c(1,-1)*threshold, linetype=3)
          + xlab(expression(paste(log[2], "(fold change)")))
          + ylab(expression(paste(-10 %.% log[10](pvalue))))
          + scale_colour_discrete(expression("Maximal difference observed by"))
          + facet_wrap(~seen, ncol=1)
          + theme(axis.text = element_text(colour = "black"))
          )
    if (.plot) {
        return(p + geom_point(aes(colour=seen)))
    } else {
        return(p)
    }
}

###########################################################
## Select for false discovery rate by Benjamini-Hochberg ##
###########################################################

reject <- function(pvec, fdr) {
    # for a vector of p-values returns the positions of tests
    # which can be rejected for the given false discovery threshold.
    order(pvec)[1:max(which(sort(pvec) <= 1:length(pvec) * fdr / length(pvec)))]
}

# accused <- function(res, fdr=0.01) {
#    return(res[reject(res$pvalue, fdr=fdr),])
#}

#' Result filtering
#' @param res result frame of test results
#' @param fdr fase discovery rate
#' @export
selectByFDR <- function (res, fdr=0.01) {
    return(droplevels(ddply(res, .(factors), function(d) d[reject(d$pvalue, fdr=fdr),])))
}

################################################
## Select by biologically feasible thresholds ##
################################################

#' Result filtering -- Test for biological effect
#'
#' The result file filtered by contains on the confidence intervals.  This
#' function will use these confidence intervals to filter out biological
#' irrelevant effects.
#' @param res Result file
#' @param threshold Biologically reasonable threshold
#' @param foldchange Is the threshold given a fold change or a log2-fold
#' change.  Default ist TRUE.
#' @export
selectByConfidence <- function (res, threshold, foldchange=TRUE) {
    if (foldchange) {
        res[(res$upr < -log2(threshold)) | (res$lwr > log2(threshold)),]
    } else {
        res[(res$upr < -threshold) | (res$lwr > threshold),]
    }
}

#' Result filtering -- Test for biological effect
#'
#' @param res Result file
#' @param cutoff the cutoff to be used in the selection
#' @export
selectByEffect <- function(res, cutoff=1) {
    return(droplevels(ddply(res, .(factors), function(d) d[abs(d$diff) > cutoff,])))
}
