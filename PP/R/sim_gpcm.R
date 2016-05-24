#' Simulate data for the gpcm model
#'
#' This function returns an integer matrix of simulated responses under given item and person parameters.
#'
#'
#'
#'@param thres An numeric matrix which contains threshold parameters for each item. The first row must contain zeroes only!
#'@param alpha A numeric vector, which contains the slope parameters - one parameter per item is expected.
#'@param theta A numeric vector which contains person parameters.
#'
#'@example ./R/.examples_simgpcm.R
#'
#'@seealso \link{sim_4pl}, \link{PP_gpcm}
#'
#'@export
#'
sim_gpcm <- function(thres, alpha, theta)
{

numbcat <- apply(thres,2,function(x) length(x)-sum(is.na(x)) - 1)  
nitem <- ncol(thres)
whereNna <- apply(thres,2,function(x) !is.na(x))

  
mat01 <- sapply(theta, function(th)
    {
      
      
    persX <- sapply(1:nitem, function(it)
        {
        thr <- thres[whereNna[,it],it]
        
        prob <- sapply(0:numbcat[it],function(ct) P_gpcm(thr,alpha[it],th,ct))
        sample(0:numbcat[it],1,prob=prob)
        
        })

    persX  
    })  

t(mat01)
}


#THETA <- rnorm(2000)


#system.time(test <- sim_gpcm(thres = THRES1,alpha = sl,theta = sort(THETA)))
# slow! :-/

