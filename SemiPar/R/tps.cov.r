########## R-function: tps.cov ##########

# Evaluates the thin plate spline
# covariance function for two dimensional
# smoothing/kriging.

# Last changed: 25 OCT 2005

tps.cov <- function(r,m=2,d=1)
{     
    r <- as.matrix(r)
    num.row <- nrow(r)
    num.col <- ncol(r)

    r <- as.vector(r)

    nzi <- (1:length(r))[r!=0]

    ans <- rep(0,length(r))

    if ((d+1)%%2!=0)    
       ans[nzi] <- (abs(r[nzi]))^(2*m-d)*log(abs(r[nzi]))     # d is even
    else
       ans[nzi] <- (abs(r[nzi]))^(2*m-d)

    if (num.col>1) ans <- matrix(ans,num.row,num.col)  # d is odd

    return(ans)
}

######### End of R-function tps.cov ########
