#' Produce a ROC curve with gini coefficient title
#'
#' This function uses ggplot to produce a themed Receiver Operator Curve and
#' calculates a Gini coefficient based on it.
#'
#' @param pred Logit/scores/probabilities to be compared against actuals
#' @param act This should be a column containing outcomes in a boolean form either as a factor or number
#' 
#' @keywords gini roc AUROC 
#' @seealso \code{AUC} \code{roc} \code{\link{giniCoef}}
#' @family creditrisk
#' @export
#' 
#' @examples 
#'   sampledata<- data.frame(val= rnorm(100) , outcome=rbinom(100,1,.8))
#'   giniChart(sampledata$val,sampledata$outcome)
#'   

giniChart <- function(pred, act) {
    stopifnot(is.numeric(pred), nlevels(factor(act)) <= 2)
    
    fpr.df <- tpr.df <- NULL  # Setting the variables to NULL first for CRAN check NOTE
    optiplum <- rgb(red = 129, green = 61, blue = 114, maxColorValue = 255)
    
    act <- factor(act)
    data <- AUC::roc(pred, act)
    
    coef <- 2 * (AUC::auc(data) - 0.5)
    
    ginidata <- data.frame(fpr.df = data$fpr, tpr.df = data$tpr)
    
    ggplot2::ggplot(ginidata, aes(x = fpr.df, y = tpr.df, colour = rgb(red = 129, green = 61, blue = 114, maxColorValue = 255))) + theme_optimum() + geom_line() + 
        scale_colour_identity() + geom_line(aes(x = fpr.df, y = fpr.df, colour = "grey")) + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + 
        labs(x = "1-Specificity", y = "Sensitivity", title = paste0("Gini = ", percent(coef)))
    
} 
