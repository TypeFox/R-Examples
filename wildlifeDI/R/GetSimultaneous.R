# ---- roxygen documentation ----
#
#' @title Identify simultaneous fixes between trajectories
#'
#' @description
#'   The function \code{GetSimultaneous} identifies and extracts simultaneous fixes, 
#'   within a given tolerance limit, between two movement datasets.
#'
#' @details
#'   This function is used to determine the simultaneous fixes between two movement 
#'   datasets facilitating further analysis.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes. For simplicity, \code{tc}
#'    is always taken in seconds.
#'
#' @return
#' A single ltraj object containing two bursts, representing the two original \code{ltraj} 
#' objects, each containing only those fixes that are deemed simultaneous.
#'
# @references
#'
#' @keywords processing
#' @seealso as.ltraj
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes
#' trajs <- GetSimultaneous(deer37, deer38, tc = 7.5*60)
#' deer37 <- trajs[1]
#' deer38 <- trajs[2]
#' 
#' @export
#
# ---- End of roxygen documentation ----
GetSimultaneous <- function(traj1,traj2,tc=0){
  
  #check to see if traj1 is not a type II
  if (attributes(traj1)$typeII == FALSE) {stop("The traj1 object is not a TypeII ltraj object")}
  #check to see if traj2 is not a type II
  if (attributes(traj2)$typeII == FALSE) {stop("The traj2 object is not a TypeII ltraj object")}
  #store as dataframes
  tr1 <- ld(traj1)
  tr2 <- ld(traj2)
  #get the length of each trajectory
  n1 <- dim(tr1)[1]
  n2 <- dim(tr2)[1]
  #Get nearest fixes that are within tc from one another.
  match1 <- data.frame()
  for (i in 1:n1){
    matched <- which.min(abs(difftime(tr1$date[i],tr2$date,units="secs")))
    temp <- data.frame(tr1=i,tr2=matched,dt=abs(difftime(tr1$date[i],tr2$date[matched],units="secs"))  )
    match1 <- rbind(match1,temp)
  }
  match2 <- data.frame()
  for (i in 1:n2){
    matched <- which.min(abs(difftime(tr2$date[i],tr1$date,units="secs")))
    temp <- data.frame(tr1=matched,tr2=i,dt=abs(difftime(tr2$date[i],tr1$date[matched],units="secs"))  )
    match2 <- rbind(match2,temp)
  }
  
  match1 <- match1[which(match1$dt <= tc),]
  match2 <- match2[which(match2$dt <= tc),]
  
  ind <- NULL
  for (i in unique(match1$tr2)){
    ind1 <- which(match1$tr2 == i)
    ind1 <- ind1[which.min(match1$dt[ind1])]
    ind <- c(ind,ind1)
  }
  match1 <- match1[ind,]
  ind <- NULL
  for (i in unique(match2$tr1)){
    ind1 <- which(match2$tr1 == i)
    ind1 <- ind1[which.min(match2$dt[ind1])]
    ind <- c(ind,ind1)
  }
  match2 <- match2[ind,]
  
  ind.1 <- which(is.na(match(match1$tr1,match2$tr1))==TRUE)
  if (length(ind.1) > 0){ match1 <- match1[-ind.1,] }
  ind.2 <- which(is.na(match(match2$tr1,match1$tr1))==TRUE)
  if (length(ind.1) > 0){ match2 <- match2[-ind.2,] }
  
  tr1.sim <- tr1[match1$tr1,]
  tr2.sim <- tr2[match1$tr2,]
  
  #convert to ltraj objects
  out.traj1 <- dl(tr1.sim)
  out.traj2 <- dl(tr2.sim)
  
  #Return the two ltraj objects
  return(c(out.traj1,out.traj2))
}
#=============== End of GetSimultaneous Function ===============================
