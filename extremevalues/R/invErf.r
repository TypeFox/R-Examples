# calculates inverse error function based on the inverse 
# cumulative normal distribution (\Phi function).
# Mark van der Loo 
#
# 20.08.2009    version 1.0
# 17.09.2009    Changed input check to allow vector input

invErf <- function(x)
{
    if ( sum(x >= 1) > 0  | sum(x <= -1) > 0 )
        stop("Argument must be between -1 and 1")

    return(qnorm((1+x)/2)/sqrt(2));

}
