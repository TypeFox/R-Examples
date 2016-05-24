#' Sensitivity analysis to correct for unknown or unmeasured polychotomous confounding without effect modification
#'
#' Simple sensitivity analysis to correct for unknown or unmeasured polychotomous
#' (3-level) confounding without effect modification. Implementation for ratio
#' measures (relative risk -- RR, or odds ratio -- OR) and difference measures
#' (risk difference -- RD).
#'
#' @param case Outcome variable. If a variable, this variable is tabulated against.
#' @param exposed Exposure variable.
#' @param implement Deprecated. Please use type instead.
#' @param type Choice of implementation, with no effect measure modification for
#' ratio measures (relative risk -- RR; odds ratio -- OR) or difference measures
#' (risk difference -- RD).
#' @param p Numeric vector defining the prevalence of the confounder. This vector
#' has 4 elements between 0 and 1, in the following order:
#' \enumerate{
#' \item Prevalence of the highest level confounder among the exposed,
#' \item Prevalence of the highest level confounder among the unexposed,
#' \item Prevalence of the mid-level confounder among the exposed, and
#' \item Prevalence of the mid-level confounder among the unexposed.
#' }
#' @param RR.cd Vector defining the confounder-disease relative risk. This vector
#' has two elements in the following order:
#' \enumerate{
#' \item Relative risk of having the highest-level confounder in diseased, and
#' \item Relative risk of having the mid-level confounder in diseased.
#' }
#' @param OR.cd Vector defining the confounder-disease odds ratio. This vector has
#' two elements in the following order:
#' \enumerate{
#' \item Odds ratio of having the highest-level confounder in diseased, and
#' \item Odds ratio of having the mid-level confounder in diseased.
#' }
#' @param RD.cd Vector defining the confounder-disease risk difference. This
#' vector has two elements in the following order:
#' \enumerate{
#' \item Risk difference of having the highest-level confounder in diseased, and
#' \item Risk difference of having the mid-level confounder in diseased.
#' }
#' @param alpha Significance level.
#' @param dec Number of decimals in the printout.
#' @param print A logical scalar. Should the results be printed?
#' 
#' @return A list with elements:
#' \item{obs.data}{The analysed 2 x 2 table from the observed data.}
#' \item{cfder1.data}{The same table for Mid-level Confounder +.}
#' \item{cfder2.data}{The same table for Highest-level Confounder +.}
#' \item{nocfder.data}{The same table for Confounder -.}
#' \item{obs.measures}{A table of relative risk with confidence intervals; Total
#' and by confounders.}
#' \item{adj.measures}{A table of Standardized Morbidity Ratio and Mantel-Haenszel
#' estimates.}
#' \item{bias.parms}{Input bias parameters.}
#'
#' @references Lash, T.L., Fox, M.P, Fink, A.K., 2009 \emph{Applying Quantitative
#' Bias Analysis to Epidemiologic Data}, pp.59--78, Springer.
#' 
#' @examples
#' # The data for this example come from:
#' # Tyndall M.W., Ronald A.R., Agoki E., Malisa W., Bwayo J.J., Ndinya-Achola J.O.
#' # et al.
#' # Increased risk of infection with human immunodeficiency virus type 1 among
#' # uncircumcised men presenting with genital ulcer disease in Kenya.
#' # Clin Infect Dis 1996;23:449-53.
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RR",
#' p = c(.6, .05, .2, .2),
#' RR.cd = c(.4, .8))
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "OR",
#' p = c(.6, .05, .2, .2),
#' OR.cd = c(.4, .8))
#' confounders.poly(matrix(c(105, 85, 527, 93),
#' dimnames = list(c("HIV+", "HIV-"), c("Circ+", "Circ-")),
#' nrow = 2, byrow = TRUE),
#' type = "RD",
#' p = c(.6, .05, .2, .2),
#' RD.cd = c(-.4, -.2))
#' @export
#' @importFrom stats qnorm
confounders.poly <- function(case,
                             exposed,
                             implement = c("RR", "OR", "RD"),
                             type = c("RR", "OR", "RD"),
                             p = NULL,
                             RR.cd = NULL,
                             OR.cd = NULL,
                             RD.cd = NULL,
                             alpha = 0.05,
                             dec = 4,
                             print = TRUE){
    if (!missing(implement)) {
        warning("Argument implement is deprecated; please use type instead.", 
                call. = FALSE)
        type <- implement
  }
    if(length(type) > 1)
        stop('Choose between RR, OR, or RD implementation.')
    if(type == "RR" & (!is.null(OR.cd) | !is.null(RD.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    if(type == "OR" & (!is.null(RR.cd) | !is.null(RD.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    if(type == "RD" & (!is.null(RR.cd) | !is.null(OR.cd)))
        stop('Mismatch between implementation type and confounder risk.')
    
    if(is.null(p))
        p <- c(0, 0)
    else p <- p
    if(length(p) != 4)
        stop('The argument p should be made of the following components: (1) Prevalence of the confounder (highest level) among the exposed, (2) Prevalence of the confounder (highest level) among the unexposed, (3) Prevalence of the confounder (mid-level) among the exposed, and (4) Prevalence of the confounder (mid-level) among the unexposed.')
    if(!all(p >= 0 & p <=1))
        stop('Prevalences should be between 0 and 1.')
    if(p[1] + p[3] >= 1)
        stop('Sum of prevalences among the exposed >= 1.')
    if(p[2] + p[4] >= 1)
        stop('Sum of prevalences among the unexposed >= 1.')
    
    if(is.null(RR.cd))
        RR.cd <- c(1, 1)
    else RR.cd <- RR.cd
    if(length(RR.cd) > 2)
        stop('Confounder-disease relative risk: more than 2 arguments.')
    if(!all(RR.cd > 0))
        stop('Confounder-disease relative risks should be greater than 0.')

    if(is.null(OR.cd))
        OR.cd <- c(1, 1)
    else OR.cd <- OR.cd
    if(length(OR.cd) > 2)
        stop('Confounder-disease odds ratio: more than 2 arguments.')
    if(!all(OR.cd > 0))
        stop('Confounder-disease odds ratios should be greater than 0.')

    if(is.null(RD.cd))
        RD.cd <- c(1, 1)
    else RD.cd <- RD.cd
    if(length(RD.cd) > 2)
        stop('Confounder-disease risk difference: more than 2 arguments.')    

    if(inherits(case, c("table", "matrix")))
        tab <- case
    else tab <- table(case, exposed)
    tab <- tab[1:2, 1:2]

    a <- tab[1, 1]
    b <- tab[1, 2]
    c <- tab[2, 1]
    d <- tab[2, 2]

    type <- match.arg(type)
    if (type == "RR") {
        crude.rr <- (a/(a + c)) / (b/(b + d))
        se.log.crude.rr <- sqrt((c/a) / (a+c) + (d/b) / (b+d))
        lci.crude.rr <- exp(log(crude.rr) - qnorm(1 - alpha/2) * se.log.crude.rr)
        uci.crude.rr <- exp(log(crude.rr) + qnorm(1 - alpha/2) * se.log.crude.rr)

        M2 <- (a + c) * p[1]
        M1 <- (a + c) * p[3]
        N2 <- (b + d) * p[2]
        N1 <- (b + d) * p[4]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- (M0 * a) / (RR.cd[2] * M1 + M0 + RR.cd[1] * M2)
        B0 <- (N0 * b) / (RR.cd[2] * N1 + N0 + RR.cd[1] * N2)
        A1 <- RR.cd[2] * A0 * M1 / M0
        B1 <- RR.cd[2] * B0 * N1 / N0
        A2 <- RR.cd[1] * A0 * M2 / M0
        B2 <- RR.cd[1] * B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRrr <- a / ((M2 * B2/N2) + (M1 * B1/N1) + (M0 * B0/N0))
        MHrr <- (A2 * N2 / (M2 + N2) + A1 * N1 / (M1 + N1) + A0 * N0 / (M0 + N0)) /
            (B2 * M2 / (M2 + N2) + B1 * M1 / (M1 + N1) + B0 * M0 / (M0 + N0))
        cfder2.rr <- (A2/(A2 + C2)) / (B2/(B2 + D2))
        cfder1.rr <- (A1/(A1 + C1)) / (B1/(B1 + D1))
        nocfder.rr <- (A0/(A0 + C0)) / (B0/(B0 + D0))
        RRadj.smr <- crude.rr / SMRrr
        RRadj.mh <- crude.rr / MHrr

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab.nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab.nocfder) <- colnames(tab)
        }
        if (print) 
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        if (print) 
            print(round(tab, dec))
        if (print)
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.rr, lci.crude.rr, uci.crude.rr))
        rownames(rmat) <- c("                       Crude Relative Risk:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.",
                                           sep = ""), "interval")
        rmatc <- rbind(c(SMRrr, RRadj.smr), c(MHrr, RRadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio", "Mantel-Haenszel")
        colnames(rmatc) <- c("SMR_RR/MH_RR", "RRc")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Relative Risk, Confounder +, Highest Level:", round(cfder2.rr, dec), "\n    Relative Risk, Confounder +, Mid-Level:", round(cfder1.rr, dec), "\n               Relative Risk, Confounder -:", round(nocfder.rr, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("Standardized Morbidity Ratio", "    SMRrr:", round(SMRrr, dec), "   RR adjusted using SMR estimate:", round(RRadj.smr, dec),
                "\nMantel-Haenszel", "                 MHrr:", round(MHrr, dec), "   RR adjusted using MH estimate:", round(RRadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("p(Confounder+HighestLevel|Exposure+):", p[1],
                "\np(Confounder+HighestLevel|Exposure-):", p[2],
                "\n    p(Confounder+MidLevel|Exposure+):", p[3],
                "\n    p(Confounder+MidLevel|Exposure-):", p[4],
                "\n  RR(ConfounderHighestLevel-Outcome):", RR.cd[1],
                "\n      RR(ConfounderMidLevel-Outcome):", RR.cd[2],
                "\n")
        rmat <- rbind(rmat, c(cfder1.rr, NA, NA), c(cfder2.rr, NA, NA),
                      c(nocfder.rr, NA, NA))
        rownames(rmat) <- c("                       Crude Relative Risk:",
                            "Relative Risk, Confounder +, Highest Level:",
                            "    Relative Risk, Confounder +, Mid-Level:",
                            "               Relative Risk, Confounder -:")
    }

    if (type == "OR"){
        crude.or <- (a/b) / (c/d)
        se.log.crude.or <- sqrt(1/a + 1/b + 1/c + 1/d)
        lci.crude.or <- exp(log(crude.or) - qnorm(1 - alpha/2) * se.log.crude.or)
        uci.crude.or <- exp(log(crude.or) + qnorm(1 - alpha/2) * se.log.crude.or)

        C2 <- c * p[1]
        C1 <- c * p[3]
        D2 <- d * p[2]
        D1 <- d * p[4]
        C0 <- c - C2 - C1
        D0 <- d - D2 - D1
        A0 <- (C0 * a) / (OR.cd[2] * C1 + OR.cd[1] * C2 + C0)
        B0 <- (D0 * b) / (OR.cd[2] * D1 + OR.cd[1] * D2 + D0)
        A1 <- OR.cd[2] * A0 * C1 / C0
        B1 <- OR.cd[2] * B0 * D1 / D0
        A2 <- OR.cd[1] * A0 * C2 / C0
        B2 <- OR.cd[1] * B0 * D2 / D0
        M2 <- A2 + C2
        N2 <- B2 + C2
        M1 <- A1 + C1
        N1 <- B1 + D1
        M0 <- A0 + C0
        N0 <- B0 + C0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        SMRor <- a / ((C2 * B2/D2) + (C1 * B1/D1) + (C0 * B0/D0))
        MHor <- (A2 * D2 / N2 + A1 * D1 / N1 + A0 * D0 / N0) /
            (B2 * C2 / N2 + B1 * C1 / N1 + B0 * C0 / N0)
        cfder2.or <- (A2 / C2) / (B2 / D2)
        cfder1.or <- (A1 / C1) / (B1 / D1)
        nocfder.or <- (A0 / C0) / (B0 / D0)
        ORadj.smr <- crude.or / SMRor
        ORadj.mh <- crude.or / MHor

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab.nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab.nocfder) <- colnames(tab)
        }
        if (print) 
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        if (print) 
            print(round(tab, dec))
        if (print)
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.or, lci.crude.or, uci.crude.or))
        rownames(rmat) <- c("                       Crude Odds Ratio:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        rmatc <- rbind(c(SMRor, ORadj.smr), c(MHor, ORadj.mh))
        rownames(rmatc) <- c("Standardized Morbidity Ratio", "Mantel-Haenszel")
        colnames(rmatc) <- c("SMR_OR/MH_OR", "ORc")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Odds Ratio, Confounder +, Highest Level:", round(cfder2.or, dec), "\n    Odds Ratio, Confounder +, Mid-Level:", round(cfder1.or, dec), "\n               Odds Ratio, Confounder -:", round(nocfder.or, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("Standardized Morbidity Ratio", "    SMRor:", round(SMRor, dec), "   OR adjusted using SMR estimate:", round(ORadj.smr, dec),
                "\nMantel-Haenszel", "                  MHor:", round(MHor, dec), "    OR adjusted using MH estimate:", round(ORadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("p(Confounder+HighestLevel|Exposure+):", p[1],
                "\np(Confounder+HighestLevel|Exposure-):", p[2],
                "\n    p(Confounder+MidLevel|Exposure+):", p[3],
                "\n    p(Confounder+MidLevel|Exposure-):", p[4],
                "\n  OR(ConfounderHighestLevel-Outcome):", OR.cd[1],
                "\n      OR(ConfounderMidLevel-Outcome):", OR.cd[2],
                "\n")
        rmat <- rbind(rmat, c(cfder1.or, NA, NA), c(cfder2.or, NA, NA),
                      c(nocfder.or, NA, NA))
        rownames(rmat) <- c("                       Crude Odds Ratio:",
                            "Odds Ratio, Confounder +, Highest Level:",
                            "    Odds Ratio, Confounder +, Mid-Level:",
                            "               Odds Ratio, Confounder -:")
    }

    if (type == "RD"){
        crude.rd <- (a / (a + c)) - (b / (b + d))
        se.log.crude.rd <- sqrt((a * c) / (a + c)^3 + (b * d) / (b + d)^3)
        lci.crude.rd <- crude.rd - qnorm(1 - alpha/2) * se.log.crude.rd
        uci.crude.rd <- crude.rd + qnorm(1 - alpha/2) * se.log.crude.rd

        M2 <- (a + c) * p[1]
        M1 <- (a + c) * p[3]
        N2 <- (b + d) * p[2]
        N1 <- (b + d) * p[4]
        M0 <- a + c - M2 - M1
        N0 <- b + d - N2 - N1
        A0 <- M0 * (a - M2 * RD.cd[1] - M1 * RD.cd[2]) / (a + c)
        B0 <- N0 * (b - N2 * RD.cd[1] - N1 * RD.cd[2]) / (b + d)
        A1 <- M1 * RD.cd[2] + A0 * M1 / M0
        B1 <- N1 * RD.cd[2] + B0 * N1 / N0
        A2 <- M2 * RD.cd[1] + A0 * M2 / M0
        B2 <- N2 * RD.cd[1] + B0 * N2 / N0
        C2 <- M2 - A2
        D2 <- N2 - B2
        C1 <- M1 - A1
        D1 <- N1 - B1
        C0 <- M0 - A0
        D0 <- N0 - B0

        if(A2 < 0 | B2 < 0 | C2 < 0 | D2 < 0 |
           A1 < 0 | B1 < 0 | C1 < 0 | D1 < 0 |
           A0 < 0 | B0 < 0 | C0 < 0 | D0 < 0)
            stop('Parameters chosen lead to negative cell(s) in adjusted 2x2 table(s).')
        
        tab.cfder2 <- matrix(c(A2, B2, C2, D2), nrow = 2, byrow = TRUE)
        tab.cfder1 <- matrix(c(A1, B1, C1, D1), nrow = 2, byrow = TRUE)
        tab.nocfder <- matrix(c(A0, B0, C0, D0), nrow = 2, byrow = TRUE)

        MHrd <- (((A2 * N2 - B2 * M2) / (M2 + N2)) +
                     ((A1 * N1 - B1 * M1) / (M1 + N1)) +
                         (((A0 * N0 - B0 * M0) / (M0 + N0)))) /
                             ((M2 * N2 / (M2 + N2)) + (M1 * N1 / (M1 + N1)) +
                                  (M0 * N0 / (M0 + N0)))
        cfder2.rd <- (A2 / M2) - (B2 / N2)
        cfder1.rd <- (A1 / M1) - (B1 / N1)
        nocfder.rd <- (A0 / M0) - (B0 / N0)
        RDadj.mh <- crude.rd - MHrd

        if (is.null(rownames(tab)))
            rownames(tab) <- paste("Row", 1:2)
        if (is.null(colnames(tab)))
            colnames(tab) <- paste("Col", 1:2)
        if (is.null(rownames(tab))){
            rownames(tab.cfder2) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder2) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder2) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder2) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.cfder1) <- paste("Row", 1:2)
        } else {
            rownames(tab.cfder1) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.cfder1) <- paste("Col", 1:2)
        } else {
            colnames(tab.cfder1) <- colnames(tab)
        }
        if (is.null(rownames(tab))){
            rownames(tab.nocfder) <- paste("Row", 1:2)
        } else {
            rownames(tab.nocfder) <- row.names(tab)
        }
        if (is.null(colnames(tab))){ 
            colnames(tab.nocfder) <- paste("Col", 1:2)
        } else {
            colnames(tab.nocfder) <- colnames(tab)
        }
        if (print) 
            cat("Observed Data:",
                "\n--------------", 
                "\nOutcome   :", rownames(tab)[1],
                "\nComparing :", colnames(tab)[1], "vs.", colnames(tab)[2], "\n\n")
        if (print) 
            print(round(tab, dec))
        if (print)
            cat("\nData, Counfounder +, Highest Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder2, dec))
        if (print)
            cat("\nData, Counfounder +, Mid-Level:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.cfder1, dec))
        if (print)
            cat("\nData, Counfounder -:",
                "\n--------------------\n\n")
        if (print)
            print(round(tab.nocfder, dec))
        if (print) 
            cat("\n")
        rmat <- rbind(c(crude.rd, lci.crude.rd, uci.crude.rd))
        rownames(rmat) <- c("                       Crude Risk Difference:")
        colnames(rmat) <- c("     ", paste(100 * (1 - alpha), "% conf.", 
                                           sep = ""), "interval")
        rmatc <- rbind(c(MHrd, RDadj.mh))
        rownames(rmatc) <- "Mantel-Haenszel"
        colnames(rmatc) <- c("MH_RD", "RDc")
        if (print)
            cat("Crude and Unmeasured Confounder Specific Measures of Exposure-Outcome Relationship:",
                "\n-----------------------------------------------------------------------------------\n\n")
        if (print) 
            print(round(rmat, dec))
        if (print)
            cat("Risk Difference, Confounder +, Highest Level:", round(cfder2.rd, dec), "\n    Risk Difference, Confounder +, Mid-Level:", round(cfder1.rd, dec),  "\n               Risk Difference, Confounder -:", round(nocfder.rd, dec), "\n")
        if (print)
            cat("\nExposure-Outcome Relationship Adjusted for Confounder:",
                "\n------------------------------------------------------\n\n")
        if (print)
            cat("\nMantel-Haenszel", "               MHrd:", round(MHrd, dec), "   RD adjusted using MH estimate:", round(RDadj.mh, dec), "\n")
        if (print)
            cat("\nBias Parameters:",
                "\n----------------\n\n")
        if (print)
            cat("p(Confounder+HighestLevel|Exposure+):", p[1],
                "\np(Confounder+HighestLevel|Exposure-):", p[2],
                "\n    p(Confounder+MidLevel|Exposure+):", p[3],
                "\n    p(Confounder+MidLevel|Exposure-):", p[4],
                "\n  RD(ConfounderHighestLevel-Outcome):", RD.cd[1],
                "\n      RD(ConfounderMidLevel-Outcome):", RD.cd[2],
                "\n")
        rmat <- rbind(rmat, c(cfder1.rd, NA, NA), c(cfder2.rd, NA, NA),
                      c(nocfder.rd, NA, NA))
        rownames(rmat) <- c("                       Crude Risk Difference:",
                            "Risk Difference, Confounder +, Highest Level:",
                            "    Risk Difference, Confounder +, Mid-Level:",
                            "               Risk Difference, Confounder -:")
    }
    invisible(list(obs.data = tab,
                   cfder1.data = tab.cfder1, cfder2.data = tab.cfder2,
                   nocfder.data = tab.nocfder,
                   obs.measures = rmat,
                   adj.measures = rmatc, 
                   bias.parms = c(p, RR.cd)))
}
