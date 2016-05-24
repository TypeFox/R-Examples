#' Simulate epidemiological record
#' 
#' Create an epidemiological record of infectors and infectees with corresponding infection and recovery times
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param lambda a numeric vector specifying the time varying intensity lambda(t) which is passed to \code{\link{simNHP.fun}} to generate the non-homogeneous Poisson process.
#' @param duration the fixed duration of infection (default is 1)
#' @param NumCases the approximate number of infected cases (default is 50)
#'
#' @return A matrix with columns "Infectee", "Infector", "InfnTime" (infection time), "RecTime" (recovery time), "DoneFlag".
#' Suitable for using with \code{\link{getLabGenealogy}}
#' 
#' @seealso \code{\link{getLabGenealogy}}
#' 
#' @importFrom NHPoisson simNHP.fun
#' 
#' @examples
#' myepirecord <- makeEpiRecord(c(1,2,3,4), duration=2, NumCases=100)
#' 
#' @export
makeEpiRecord <- function(lambda,duration=1,NumCases=50) {
  lambda <- c(lambda,0) # required because of < instead of <= in simNHP.fun
  epirecord <- matrix(0,1,5) # infectee, infector, infection time, recovery time, done flag
  cis <- 1 # cis for current and pending infectors
  Infectee <- 1; 
  rec <- duration; # note fixed duration of infection. 
  epirecord[1,] <- c(1,0,0, rec, 0); 
  while (nrow(epirecord)<=NumCases & length(cis)>0) {
    itimes <- simNHP.fun(lambda)$posNH
    itimes <- itimes/(length(lambda)-1); # scaled so in [0,1]
    ci <- cis[1] # current infector. 
    NumInfected <- length(itimes) 
    if (NumInfected>0) {
      for (n in 1:NumInfected){
        Infectee <- Infectee+1; 
        epirecord <- rbind(epirecord,c(Infectee,ci, epirecord[ci,3]+itimes[n],epirecord[ci,3]+itimes[n]+rec,0))
        cis <- c(cis,Infectee)
      }
    }
    cis <- cis[-1] # 
  }
  colnames(epirecord) <- c("Infectee","Infector","InfnTime","RecTime","DoneFlag")
  return(epirecord)
}
