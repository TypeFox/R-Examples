#' Sensitivity analysis to correct for selection bias.
#'
#' Simple sensitivity analysis to correct for selection bias using estimates of
#' the selection proportions.
#'
#' @param case Outcome variable. If a variable, this variable is tabulated
#' against.
#' @param exposed Exposure variable.
#' @param selprob Numeric vector defining the selection probabilities. This vector has 4 elements between 0 and 1, in the following order:
#' \enumerate{
#' \item Selection probability among cases exposed,
#' \item Selection probability among cases unexposed,
#' \item Selection probabillity among noncases exposed, and
#' \item Selection probability among noncases unexposed.
#' }
#' @param alpha Significance level.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#' 
#' @return A list with elements:
#' \item{obs.data}{The analysed 2 x 2 table from the observed data.}
#' \item{corr.data}{The same table corrected for  selection proportions.}
#' \item{obs.measures}{A table of odds ratios and relative risk with confidence intervals.}
#' \item{adj.measures}{Selection bias corrected measures of outcome-exposure relationship.}
#' \item{bias.parms}{Input bias parameters: selection probabilities.}
#'
#' @examples
#' # The data for this example come from:
#' # Stang A., Schmidt-Pokrzywniak A., Lehnert M., Parkin D.M., Ferlay J., Bornfeld N.
#' # et al.
#' # Population-based incidence estimates of uveal melanoma in Germany. Supplementing
#' # cancer registry data by case-control data.
#' # Eur J Cancer Prev 2006;15:165-70.
#' selection(matrix(c(136, 107, 297, 165),
#' dimnames = list(c("UM+", "UM-"), c("Mobile+", "Mobile-")),
#' nrow = 2, byrow = TRUE),
#' selprob = c(.94, .85, .64, .25))
#' @export
#' @importFrom stats qnorm
selection <- function(case,
                      exposed,
                      selprob = NULL,
                      alpha = 0.05,
                      dec = 4,
                      print = TRUE) {
    if(is.null(selprob))
        selprob <- c(1, 1, 1, 1)
    else selprob <- selprob
    if(!is.vector(selprob))
        stop('The argument selprob should be a vector of length 4')
    if(length(selprob) != 4)
        stop('The argument selprob should be made of 4 components in the following order: (1) Selection probability among cases exposed, (2) Selection probability among cases unexposed, (3) Selection probability among noncases exposed, and (4) Selection probability among noncases unexposed')
    if(!all(selprob >= 0 & selprob <=1))
        stop('Selection probabilities should be between 0 and 1')

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else tab <- table(case, exposed)
    tab <- tab[1:2, 1:2]

    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    rr <- (a/(a + c)) / (b/(b + d))
    se.log.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
    lci.rr <- exp(log(rr) - qnorm(1 - alpha/2) * se.log.rr)
    uci.rr <- exp(log(rr) + qnorm(1 - alpha/2) * se.log.rr)

    or <- (a/b) / (c/d)
    se.log.or <- sqrt(1/a + 1/b + 1/c + 1/d)
    lci.or <- exp(log(or) - qnorm(1 - alpha/2) * se.log.or)
    uci.or <- exp(log(or) + qnorm(1 - alpha/2) * se.log.or)

    A0 <- a / selprob[1]
    B0 <- b / selprob[2]
    C0 <- c / selprob[3]
    D0 <- d / selprob[4]

    tab.corr <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)
    rr.corr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
    or.corr <- (A0/B0) / (C0/D0)
   
    if (is.null(rownames(tab)))
        rownames(tab) <- paste("Row", 1:2)
    if (is.null(colnames(tab)))
        colnames(tab) <- paste("Col", 1:2)
    if (is.null(rownames(tab))){
        rownames(tab.corr) <- paste("Row", 1:2)
        } else {
        rownames(tab.corr) <- row.names(tab)
    }
    if (is.null(colnames(tab))){ 
        colnames(tab.corr) <- paste("Col", 1:2)
        } else {
        colnames(tab.corr) <- colnames(tab)
    }
    if (print) 
        cat("Observed Data:", "\n---------------------------------------------------", 
            "\nOutcome   :", rownames(tab)[1],
            "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
    if (print) 
        print(round(tab, dec))
    if (print)
        cat("\nData Corrected for Selected Proportions:", "\n---------------------------------------------------\n\n")
    if (print)
        print(round(tab.corr, dec))
    if (print) 
        cat("\n")
    rmat <- rbind(c(rr, lci.rr, uci.rr), c(or, lci.or, uci.or))
    rownames(rmat) <- c("Observed Relative Risk:", "   Observed Odds Ratio:")
    colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
        sep = ""), "interval")
    rmatc <- rbind(rr.corr, or.corr)
    rownames(rmatc) <- c("Selection Bias Corrected Relative Risk:",
                         "   Selection Bias Corrected Odds Ratio:")
    if (print) 
        print(round(rmat, dec))
    if (print)
        cat("\n")
    if (print)
        print(round(rmatc, dec))
    if (print)
        cat("\n")
    bias <- rbind(selprob[1], selprob[2], selprob[3], selprob[4])
    rownames(bias) <- c("     Selection probability among cases exposed:",
                        "   Selection probability among cases unexposed:",
                        "  Selection probability among noncases exposed:",
                        "Selection probability among noncases unexposed:")
    if (print)
        print(bias)
    invisible(list(obs.data = tab, corr.data = tab.corr,
                   obs.measures = rmat, adj.measures = rmatc,
                   bias.parms = bias))
}
