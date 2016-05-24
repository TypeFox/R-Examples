test.deleteViews <- function()
{    
    stocks <- colnames(monthlyReturns)
    viewInput <- matrix(0, ncol = 6, nrow = 2, dimnames = list(NULL, c(stocks)))

    viewInput[1,"IBM"] <- 1
    viewInput[1, "DELL"] <- -1
    viewInput[2, "C"] <- 1
    viewInput[2, "JPM"] <- -1
    confidences <- c(70, 10) 
    
    x <- BLViews(viewInput, c(-0.1, 0.1), confidences, stocks)
    x <- deleteViews(x,1)
    y <- BLViews(viewInput[-1,,drop=FALSE], c(0.1), confidences[2], stocks)
    checkEquals(x, y)
    
    viewInput <- matrix(0, ncol = 6, nrow = 2, dimnames = list(NULL,stocks))
    viewInput[1,"IBM"] <- 1
    viewInput[1, "DELL"] <- -1
    viewInput[2, "C"] <- 1
    viewInput[2, "JPM"] <- -1
    confidences <- c(0.7, 0.1)
    viewDist <- list(distribution("unif", min = 0.04, max = 0.1), distribution("norm", sd = 10, mean = 0.05))
    
    x <- COPViews(viewInput, viewDist, confidences, stocks)
    x <- deleteViews(x,1)
    y <- COPViews(viewInput[2,,drop=FALSE], viewDist[2], confidences[2], stocks)
    checkEquals(x,y)
    
}