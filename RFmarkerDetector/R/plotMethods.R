#' Plotting the average AUC as a function of the number of combinations
#' 
#' This function allows to plot the average AUC as a function of the number k-combinations of the n input
#' variables. If n is the number of input variables, the number of k-combinations of those variables is equal to
#' \eqn{n!/k!(n!-k!)}. Each of these combinations contains the indexes of the input variables selected. For each 
#' combination we can extract a dataset, build a random forest model and perform a cross-validation. We can describe
#' the performance of each cross-validated model with an 'average' ROC curve and its AUC. The collected auc values
#' for each combination (dataset) are used by the function to build a diagram of the AUC as a function 
#' of the number of combinations
#' @param auc_values an array with the auc values to plot
#' @param num_of_variables the k dimension of each combination
#' @param num_of_combinations the number of k combinations of the set of the input variables
#' @importFrom graphics mtext
#' @author Piergiorgio Palla

plotAUCvsCombinations <- function(auc_values, num_of_variables, num_of_combinations) {
    
    x_lab = "combinations"
    y_lab = "Average AUC"
    
    plot(auc_values, type = "l", col = "green", xlab = x_lab, ylab = y_lab, ylim = c(0.3, 1))
    xmin = which.min(auc_values)
    ymin = min(auc_values)
    xmax = which.max(auc_values)
    ymax = max(auc_values)
    points(xmax, ymax, col = "red")
    abline(v = xmax, h = ymax, col = "red", lty = 3)
    text(x = (xmax + 0.05), y = (ymax + 0.01), labels = paste("AUC max: ", ymax))
    mtext(text = paste("# var:", num_of_variables), side = 3, line = 2)
    mtext(text = paste("# comb:", num_of_combinations), side = 3, line = 1)
    grid()
    
    
} 
