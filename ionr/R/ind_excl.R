#' Wrapper and scripts for indicator exclusion procedure
#'
#' @description  For each item, correlation between the scale's sum scores and outcome is
#' calculated such that the particular item is excluded from the sum scores. Each of the
#' obtained correlations will then be compared with the original scale–outcome correlation
#' (sum score of all items). This comparison can be conducted with William's test for two
#' dependent correlations that share one variable (Steiger, 1980), using
#' \code{\link[psych]{r.test}} from \pkg{psych}. William's test characterises difference
#' between correlations with a p-value—a small p-value indicates that the tested
#' difference between correlations is unlikely to have happened by chance and could be
#' considered a real difference. Thus, each item will receive a p-value characterising
#' the 'significance' of difference between correlations—here called 'Significance Of
#' iNdicator Exclusion' (SONE). When one item is excluded, another round is begun
#' until no more items should be excluded.

# simple wrapper for running the indicator exclusion procedure and plotting the results
#' @param multi influences cex of certain plot variables. Defaults to 1
#' @param indicators Set of numeric indicators (items) in a matrix.
#' @param indicators2 An additional set of indicators (e.g. informant-report )
#' @param outcome A numeric outcome vector. Indicators and outcome can be simulated with
#' \code{\link{scale_sim}}
#' @param covars A data frame with covariates to take into account. The outcome is residualised for
#' these covariates. For instance, BMI was residualised for age, gender, and education in Vainik
#' et al. (2015) EJP.
#' @param scalename A string for labelling the scale
#' @param outcomename A string for labelling the outcome
#' @param indicatornames An array of strings for labelling the outcome. Default to numbers
#'  from 1 to n of indicators

#' @param pcrit a p-value characterising the ‘significance’ of difference between
#' correlations—here called ‘significance of indicator exclusion’ (SONE). Look it up
#' from Table 2 in Vainik, Mõttus et al 2015, or simulate using
#'  \code{\link{optimal_p}} function
#' @param location1 Location for legends at left-side plot
#' @param location2 Location for legends at right-side plot

#' @param draw TRUE plots the result to a .tiff file in the working directory.
#' Defaults to FALSE
#' @param subset Allows exluding certain indicators from the start. Use numbers
#' @param coruse argument for function cor(). Defaults to 'everything', as
#' simulations have no missing data.
#' @param verbose option for observing steps for debugging. Defaults to FALSE
#' @param ci should output object and plot have 95% confidence intervals (CI-s). Defaults to
#' CI-s from  \code{\link[psych]{corr.test}}. If you insert a number (e.g., ci=5000), then the CI-s
#' are bootstrapped using \code{\link[psych]{cor.ci}}. Any other string results in no CI-s.
#' r value in output matrix is taken from \code{\link[stats]{cor}}.
#'
#' @importFrom stats as.formula cor lm quantile resid rnorm runif sd
#' @importFrom graphics axis legend lines mtext par plot segments text title
#' @importFrom utils flush.console
#' @importFrom grDevices dev.off tiff

#' @export
#'
#' @return Plots results using , using \code{\link[gplots]{barplot2}} from  \pkg{gplots}.
#' Also returns scale-outcome correlation magnitude(s) and their comparison, if
#'  appropriate
#' @encoding utf-8
#' @examples
#' ### Create a scale-outcome set that violates IOn_ Only 2 indicators out of 8
#' ### relate to the outcome, the others just relate to the 2 indicators. This setting is
#' ### similar to the N5: Impulsiveness - BMI association in Vainik et al (2015) EJP paper.
#' set.seed(466)
#' a<-scale_sim(n=2500, to_n=2, tn_n=6)
#' # Last 2 indicators have considerably higher correlation with the outcome
#' ind_excl(a[[1]], outcome=a[[2]], pcrit=0.0037)
#'
#' ## boostrapped confidence intervals
#' ind_excl(a[[1]], outcome=a[[2]], pcrit=0.0037, ci=100)
#'
#' # no confidence intervals
#' ind_excl(a[[1]], outcome=a[[2]], pcrit=0.0037, ci="no")
#'
#' ## include covariates in the model
#' covx=rnorm(2500)
#' covy=rnorm(2500)
#' outcome=a[[2]]+0.3*covx+0.4*covy
#' covars=data.frame(covx=covx, covy=covy)
#'
#' # ind_excl() with covariates taken into account
#' ind_excl(a[[1]],outcome=outcome,covars=covars, pcrit=0.0037)
#'
#' #effect sizes are lower when noisy covariatse are not accounted for
#' ind_excl(a[[1]],outcome=outcome, pcrit=0.0037)
#'
#' # just a single covariate also needs to be in data frame
#'
#' covx=rnorm(2500)
#' outcome=a[[2]]+0.3*covx
#' covars=data.frame(covx=covx)
#'
#' ind_excl(a[[1]],outcome=outcome,covars=covars, pcrit=0.0037)
#'
#' ### Create a scale-outcome set that has ION, all 8 indicators relate to the outcome
#' set.seed(466)
#' b<-scale_sim(n=2500, to_n=8, cor_to_outcome = 0.35)
#' # All indicators correlate largely on the same level with the outcome.
#' ind_excl(b[[1]], outcome=b[[2]], pcrit=1.7*10^-4)
#'
#' #note that using cor_to_outcome=0.25, sometimes still indicators get wrongly flagged.
#' # Here, the method could probably be improved..
#'
#' ### Create a scale-outcome set that violates ION - only 1 indicator relates to the
#' ### outcome. Include other-report.
#' set.seed(466)
#' c<-scale_sim(n=2500, to_n=1, tn_n=7, indicators2=TRUE)
#' # Last indicator has considerably higher correlation with the outcome
#' ind_excl(c[[1]], c[[3]], outcome=c[[2]], pcrit=0.0037)


# for debugging
# a<-scale_sim(n=2500, to_n=2, tn_n=6) indicators=a[[1]] outcome=a[[2]] indicators2 = vector() scalename =
# 'scale'; outcomename = 'outcome'; indicatornames = 1:ncol(indicators); tagged = vector(); tagged2 = vector();
# location1 = 'topleft'; location2 = 'topright'; pcrit = 0.0037; multi = 1; coruse = 'everything' subset =
# vector() draw=F


ind_excl <- function(indicators, indicators2 = vector(), outcome, covars= NULL, scalename = "scale", outcomename = "outcome", indicatornames = 1:ncol(indicators),
    pcrit, location1 = "topleft", location2 = "topright", draw = F, subset = vector(), coruse = "everything", multi = 1,
    verbose = F,ci="estimate") {

    # residualise outcome for the covariates if needed

    if(!is.null(covars)){
        f = as.formula(paste("outcome", paste(names(covars), collapse=" + "), sep=" ~ "))
        outcome = resid(lm(f, data=covars, na.action="na.exclude"))
        cat("You included co-variates: ", paste(names(covars), collapse=", "), "\n")
    }



    # exclude certain indicators from start, if wished ( by number, not by name)
    if (length(subset) != 0) {
        indicators <- indicators[, subset]
        if (length(indicators2) != 0)
            indicators2 <- indicators2[, subset]
        indicatornames <- indicatornames[subset]
    }
    # tag indicators that should be excluded
    tagged <- ind_excl_inc(indicators, outcome, indicatornames, pcrit, coruse = coruse, verbose = verbose)
    if (length(indicators2) != 0)
        tagged2 <- ind_excl_inc(indicators2, outcome, indicatornames, pcrit, coruse = coruse, verbose = verbose)

    # in case you want to plot results to file
    if (draw == T) {
        scalename2 <- gsub(": ", "-", scalename)
        tiff(paste0(scalename2, ".tiff"), res = 300, width = 30, height = 20, units = "cm", compression = "lzw", pointsize = 15)
        par(mfrow = c(1, 1), mar = c(5, 6, 4, 2))
    }

    # plot either one or two set of indicators

    if (length(indicators2) != 0) {
        # plot both

        cors=ind_excl_plot(indicators = indicators, indicators2 = indicators2, outcome = outcome, scalename = scalename,
            outcomename = outcomename, indicatornames = indicatornames, tagged = tagged, tagged2 = tagged2, location1 = location1,
            location2 = location2, pcrit = pcrit, coruse = coruse, multi = multi,ci=ci)
    } else {

        # plot just one
        cors=ind_excl_plot(indicators = indicators, indicators2 = vector(), outcome = outcome, scalename = scalename, outcomename = outcomename,
            indicatornames = indicatornames, tagged = tagged, tagged2 = vector(), location1 = location1, location2 = location2,
            pcrit = pcrit, coruse = coruse, multi = multi,ci=ci)
    }

    # finish plotting to a file
    if (draw == T)
        dev.off()
    return(cors)

}
#

