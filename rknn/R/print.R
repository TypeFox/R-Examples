################################################################################
# Random KNN Print Functions                                                   #
# File:   print.R                                                              #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
# Change Log:                                                                  #
#          March 03, 2010 -- add ... parameter                                 #
################################################################################
print.rknnBE <- function (x, ...) 
{ 
    if (!inherits(x, "rknnBE")) stop(deparse(substitute(x)), " is not a rknnBE object");
    print(cbind(p = x$p, mean_accuracy = x$mean_accuracy))
}

print.rknn <- function (x, ...) 
{
    if (!inherits(x, "rknn")) 
        stop(deparse(substitute(x)), " is not a rknn object");
    cat("\nCall:\n", strwrap(deparse(x$call), width = getOption("width")), 
        "\n", fill = TRUE, sep = "");
    cat("Number of neighbors: ", x$k, "\n");
    cat("Number of knns: ", x$r, "\n");
    cat("No. of variables used for each knn: ", x$mtry, "\n");
    cat("Prediction: \n");
    
    print.default(format(x$pred, digits = 3), print.gap = 1, quote = FALSE)
        
    invisible(x);
}
print.rknnSupport <-function (x, ...) 
{
    if (!inherits(x, "rknnSupport")) 
        stop(deparse(substitute(x)), " is not a rknnSupport object")
        
    cat("\nCall:\n", strwrap(deparse(x$call), width = getOption("width")), 
        "\n", fill = TRUE, sep = "");
    cat("Number of knns: ", x$r, "\n")
    cat("No. of variables used for each knn: ", x$mtry, "\n");
    cat("Accuracy:", x$accuracy, "\n");
    if (ncol(x$confusion) < 15) {
        cat("Confusion matrix:\n")
        print.default(format(x$confusion, digits = 3), print.gap = 1, quote = FALSE)
    };
    
    invisible(x);
}

################################################################################
