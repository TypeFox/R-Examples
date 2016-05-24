"drmLOFls" <- function()
{
    ## Defining lack-of-fit/goodness-of-fit tests
    anovaTest <- function(formula, ds)
    {
        anovaFit <- lm(formula, data = ds)
        if (df.residual(anovaFit) > 0)
        {
            return(list(test = "F", anovaFit = anovaFit))
        } else {
            return(NULL)
        }
    }

    gofTest <- NULL

    return(list(anovaTest = anovaTest, gofTest = gofTest))
}