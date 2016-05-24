"drmLOFbinomial" <- function()
{
    ## Defining a goodness-of-fit test
    gofTest <- function(resp, weights, fitted, dfres)
    {
        ## Removing 0s and 1s in fitted values
        zeroTol <- 1e-12  # no global constant
        indVec <- ( (fitted < zeroTol) | (fitted > 1-zeroTol) )
        dfReduc <- sum(indVec)

        total <- weights  # (object$"data")[, 5]
        success <- resp*weights  # total*(object$"data")[, 2]
        expected <- total*fitted  # fitted(object)

        ## Pearson's statistic (sum of squared Pearson residuals)
        c( sum( ((success - expected)^2 / (expected*(1 - fitted)))[!indVec] ), dfres - dfReduc)  # df.residual(object))
    }


    ## Defining goodness-of-fit function
    anovaTest <- function(formula, ds)
    {
#       count <- resp*weights
        anovaFit <- glm(formula, family=binomial(link = "logit"), data=ds)
        if (df.residual(anovaFit)>0)
        {
            return(list(test = "lr", anovaFit = anovaFit))
        } else {
            return(NULL)
        }
    }
    anovaTest <- NULL  # lack-of-fit test not meaningful in most situations

    return(list(anovaTest = anovaTest, gofTest = gofTest))
}
