test.addBLViews <- function()
{  
    # very simple case first
    stocks <- colnames(monthlyReturns)

    viewInput <- matrix(0, ncol = 6, nrow = 2, dimnames = list(NULL, stocks))
    viewInput[1,"IBM"] <- 1
    viewInput[1, "DELL"] <- 0.04
    viewInput[2, "C"] <- 1
    viewInput[2, "JPM"] <- 0.6

    confidences <- 1 / c(0.7, 0.1)

    views <- new("BLViews", P = viewInput, qv = c(1,1) , 
                    confidences = confidences, assets = stocks)
    x <- BLViews(P = viewInput, confidences = confidences, q = c(1,1), assetNames = stocks)
    checkEquals(x, views)

    pick <- matrix(0, ncol = 2, nrow = 1, dimnames = list(NULL, c("MS", "BAC")))
    pick[1, "MS"] <- 1
    pick[1, "BAC"] <- 0.01
        
    x <- addBLViews(pick, 0.15, 1/ 0.03 , views = x)
   
    views2 <- BLViews(P = matrix(0,nrow = 3, ncol = 6, dimnames=list(NULL,stocks)), 
                    q = rep(0, 3), confidences = 1 / c(0.7,0.1,0.03),stocks)
    
    views2@P[1:2,] <- viewInput
   
    views2@P[3, c("MS", "BAC")] <- c(1, 0.01)
  
    views2@qv <- c(1,1,0.15)
    checkEquals(x, views2)
    
}

test.addCOPViews <- function()
{

    stocks <- colnames(monthlyReturns)

    viewInput <- matrix(0, ncol = 6, nrow = 2, dimnames = list(NULL, stocks))
    viewInput[1,"IBM"] <- 1
    viewInput[1, "DELL"] <- 0.04
    viewInput[2, "C"] <- 1
    viewInput[2, "JPM"] <- 0.6
    
    confidences <- c(0.7, 0.1)
    viewDist <- list(distribution("unif", min = 0.04, max = 0.1), distribution("norm", sd = 10, mean = 0.05))
    views1 <- COPViews(viewInput, viewDist, confidences, stocks)
    views2 <- new("COPViews", pick = viewInput, viewDist = viewDist, confidences = confidences, assets = stocks)
    
    checkEquals(views1, views2)
    
    pick <- matrix(0, ncol = 2, nrow = 1, dimnames = list(NULL, c("MS", "BAC")))
    pick[1, "MS"] <- 1
    pick[1, "BAC"] <- 0.01
    
	if(!require("sn", quiet = TRUE))
	{
		warning("The next tests require the sn package, which could not be loaded \n")
		return()
	}
    viewDist2 <- list(distribution("sn", xi = 0.05, omega = 10, alpha = 0.001))
    views1 <- addCOPViews(pick, viewDist2, 0.4, views1 )
    
    views2 <- COPViews(pickMatrix = matrix(0, ncol = 6, nrow = 3, dimnames = list(NULL, stocks)), 
        c(viewDist, viewDist2), c(confidences, 0.4), stocks )
    
    views2@pick[1:2,] <- viewInput
    views2@pick[3, c("MS", "BAC")] <- c(1, 0.01)
  
    checkEquals(views1, views2)

}

test.newPMatrix <- function()
{
    x <- newPMatrix(c("DAX", "FTSE", "CAC"), 2)
    checkEquals(x, matrix(0, ncol = 3, nrow = 2, dimnames = list(NULL, c("DAX", "FTSE", "CAC"))))
}