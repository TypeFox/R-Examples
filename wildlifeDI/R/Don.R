# ---- roxygen documentation ----
#
#' @title Doncaster's measure of dynamic interaction
#'
#' @description
#' The function \code{Don} measures the dynamic interaction between two moving objects following
#' the methods outlined by Doncaster (1990).
#'
#' @details
#' This function can be used to compute the Doncaster (1990) methods for measuring
#' dynamic interaction between two objects. The Doncaster method tests the proportion
#' of simultaneous fixes that are below \code{dc} against that which would be
#' expected based on the distribution of distances between all fixes.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'         two fixes are spatially together.
#'
#' @return
#' This function first returns a plot, for distance values ranging from 0 to the maximum
#' distance separating two fixes, of the observed proportion of simultaneous fixes below
#' each distance value. The expected values based on all fixes are also included.
#' Second, a list is returned that contains the contingency table of simultaneous fixes
#' (paired) and non-paired fixes below and above \code{dc}, along with the
#' associated \emph{p}-value from the Chi-squared test.
#' \itemize{
#'  \item conTable -- contingency table showing frequency of paired and non-paired fixes above and below \code{dc}.
#'  \item p.value -- \emph{p}-value from the Chi-squared test of \code{conTable}.
#'  }
#'
#' @references
#' Doncaster, C.P. (1992) Non-parametric estimates of interaction from radio-tracking
#' data. \emph{Journal of Theoretical Biology}, \bold{143}: 431-443.
#'
#' @keywords indices
#' @seealso GetSimultaneous
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' Don(deer37, deer38, tc = 7.5*60, dc = 50)
#' @export
#
# ---- End of roxygen documentation ----
Don <- function(traj1,traj2,tc=0,dc=50){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  
  n <- nrow(tr1)

  euc <- function(x1,y1,x2,y2){sqrt(((x1 - x2)^2) + ((y1 - y2)^2))}
  
  #calculate the observed (Paired) distances
  Do <- euc(tr1$x,tr1$y,tr2$x,tr2$y)
  
  #calculate the expected distances overall
  De <- matrix(nrow=n,ncol=n)
  for (i in 1:n){
    De[i,] <- euc(tr1$x[i],tr1$y[i],tr2$x,tr2$y)
    De[i,i] <- NA
  }
  
  #compute the cumulactive frequency plot
  pd <- seq(0,max(c(Do,De),na.rm=T),length.out=50)
  Co <- rep(0,length(pd))
  Ce <- Co
  for (i in 1:length(pd)){
    Co[i] <- length(which(Do < pd[i]))
    Ce[i] <- length(which(De < pd[i])) + Co[i]
  }
  Co <- Co / n
  Ce <- Ce / (n^2)
  
  plot(pd,Ce,type="l",col="grey50",ylab="probability",xlab="Distance <= (m)")
  points(pd,Co,pch=20)
  
  #compute the contingency table and Chi-Squared test.
  con <- matrix(nrow=3,ncol=3,dimnames=list(c("Paired","Non-Paired","Totals"),c("below.crit","above.crit","Totals")))
  
  con[1,1] <- length(which(Do < dc))
  con[2,1] <- length(which(De < dc))
  con[3,1] <- length(which(Do < dc)) +  length(which(De <= dc))
  
  con[1,2] <- length(which(Do > dc))
  con[2,2] <- length(which(De > dc))
  con[3,2] <- length(which(Do > dc)) + length(which(De > dc))
  
  con[,3] <- c(n,(n^2)-n, n^2)
  
  chi.p <- chisq.test(con[1:2,1:2])$p.value
  
  #print(paste("A critical distance of ",dc," was used, resulting in a p-value = ", chi.p,".",sep=""))
  return(list(conTable=con,p.value=chi.p))
}
#======================== End of Doncaster Function ============================