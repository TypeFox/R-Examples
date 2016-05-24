#' Calculate Harrell's c-index
#' 
#' This function calculates Harrell's c-index.
#' 
#' 
#' @param formula Formula for prediction model to be used as in
#' \code{\link[survival:coxph]{coxph}}
#' @param data Data set in which to interpret the formula
#' @return A list with elements \item{concordant}{The number of concordant
#' pairs} \item{total}{The total number of pairs that can be evaluated}
#' \item{cindex}{Harrell's c-index}
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Harrell FE, Lee KL & Mark DB (1996), Multivariable prognostic
#' models: issues in developing models, evaluating assumptions and adequacy,
#' and measuring and reducing errors, Statistics in Medicine 15, 361-387.
#' 
#' van Houwelingen HC, Putter H (2012). Dynamic Prediction in Clinical Survival
#' Analysis. Chapman & Hall.
#' @keywords univar
#' @examples
#' 
#' data(ova)
#' cindex(Surv(tyears, d) ~ Karn + Broders + FIGO + Ascites + Diam, data = ova)
#' 
#' @export cindex
cindex <- function(formula, data)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    n <- length(time)
    ord <- order(time,-status)
    time <- time[ord]
    status <- status[ord]
    x <- x[ord]
    # pairs (i,j) for which the smallest observed time is an event time
    wh <- which(status==1)
    total <- concordant <- 0
    for (i in wh) {
        for (j in ((i+1):n)) {
            if (time[j] > time[i]) {# ties not counted
                total <- total + 1
                if (x[j] < x[i]) concordant <- concordant + 1
                if (x[j] == x[i]) concordant <- concordant + 0.5
            }
        }
    }
    return(list(concordant=concordant,total=total,cindex=concordant/total))
}
