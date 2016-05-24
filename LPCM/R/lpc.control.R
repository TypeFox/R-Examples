lpc.control <- function (iter =100,  cross=TRUE, boundary = 0.005, convergence.at =0.00001, mult=NULL,  ms.h=NULL, ms.sub=30,  pruning.thresh=0.0, rho0=0.4) 
{
    if (boundary!=0 &&  boundary < convergence.at) {
        warning("The boundary correction will not have any effect if its threshold is set to a smaller value than the convergence criterion.") 
      }
    if (iter <= 9) {
        warning("Please choose iter=10 at least. The default value iter=100 has been used instead.")
        iter <- 100
      }
    if (pruning.thresh<0) {
        warning("This should be a non-negative number. The default 0 has been used instead.") 
    pruning.thresh<-0
    }
    if (rho0<0) {
        warning("This should be a non-negative number. The default 0.4 has been used instead.") 
    rho0<-0.4
    }    
    if (!is.logical(cross)){
       warning("cross needs to be a Boolean. The default FALSE has been used instead.")
     }
    if (!is.null(mult)){
       if (mult<1){
          mult <-1
          warning("mult needs to be a positive integer. The value mult=1 has been used instead.")
       }    
       cat("Number of starting points manually enforced to be", mult, ".", "\n")
     }
    if ((ms.sub<0) || (ms.sub >100)){
       ms.sub<- 30
       warning("ms.sub needs to be a percentage between 0 and 100. The default choice 30 has been used instead.")
     }
    list(iter = iter, cross=cross,
         boundary = boundary, convergence.at = convergence.at,
         mult=mult, ms.h=ms.h, ms.sub=ms.sub,
         pruning.thresh=pruning.thresh, rho0=rho0)
}
