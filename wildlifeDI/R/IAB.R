# ---- roxygen documentation ----
#
#' @title Benhamou's IAB Index
#' 
#' @description
#' The function \code{IAB} computes the IAB index following the methods described in the paper by 
#' Benhamou et al. (2014). It facilitates global analysis, with the significance testing procedure
#' described in the paper, but also a local level output, to explore the IAB statistic through time.
#'
#' @details
#' The function \code{IAB} can be used to test for direct interaction in wildlife telemetry data and affords
#' a novel significance testing procedure that takes into account the serially correlated structure of
#' telemetry data. Specifically, it computes an index analogous to the Bhattacharyya coefficient between the 
#' potential influence domains of two animals. Like the other indices, IAB is dependent on the selection of an
#' appropriate value for \code{dc} (which is termed \eqn{\Delta}{Delta} in the article). The \code{dc} parameter here is 
#' not a threshold distance, but rather the distance at which the function shows maximum slope (see 
#' Benhamou et al. 2014). 
#' 
#' The significance testing procuedure uses a wrapped shifting method in order to maintain the serially correlated 
#' structure of the data. At each shift, a sample value of IAB (termed MAB) is computed in order to generate 
#' a distribution of values to test against (for more information see Benhamou et al. 2014). 
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc critical distance where the IAB function will show maximium slope -- see Benhamou et al. (2014) 
#'    for more advice on selecting this parameter.
#' @param local logical value indicating whether a dataframe (\code{local = TRUE}) containing the IAB 
#'    index for each simultaneous fix should be returned, or (if \code{local = FALSE} - the default) 
#'    the global index along with associated p-value.
#'
#' @return
#' If \code{local=FALSE} (the default) IAB returns the numeric value of the IAB index and the associated p-value.
#' If \code{local=TRUE} IAB returns a dataframe (containing the date/times of \emph{all} simultaneous fixes, 
#' along with the distance between fixes at each time, and the IAB index value for each simultaneous fix.
#'
#' @references
#' Benhamou, S., Valeix, M., Chamaille-Jammes, S., Macdonald, D., Loveridge, A.J. (2014) Movement-based analysis
#' of interactions in African lions. \emph{Animal Behaviour}, \bold{90}: 171-180.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI, Prox
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' IAB(deer37, deer38, tc=7.5*60, dc=50)
#' df <- IAB(deer37, deer38, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
# ---- End of roxygen documentation ----


### Simon Benhamou's IAB method for ST interaction
IAB <- function(traj1,traj2,tc=0,dc=50,local=FALSE){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)
  #Calculate the observed distances
  IAB.df <- data.frame(date=tr1$date,Dab=sqrt((tr1$x - tr2$x)^2 + (tr1$y - tr2$y)^2))
  IAB.df$Iab <- exp((-1/2)*(IAB.df$Dab/dc)^2)
  
  if (local == FALSE){
    IAB. <- mean(IAB.df$Iab,na.rm=TRUE)
    
    #Significance Testing
    #-------------------------
    MAB <- NULL
    for (k in 1:(n-1)){
      tr2. <- rbind(tr2[(k+1):n,],tr2[1:k,])
      MAB.df <- data.frame(Dab=sqrt((tr1$x - tr2.$x)^2 + (tr1$y - tr2.$y)^2))
      MAB.df$Mab <- exp((-1/2)*(MAB.df$Dab/dc)^2)
      MAB <- c(MAB,mean(MAB.df$Mab,na.rm=TRUE))
    }
    ng <- length(which(MAB > IAB.))
    nb <- length(which(MAB < IAB.))
    P.attract <- (ng + 1)/n
    P.avoid <- (nb + 1)/n
    #-------------------------
    
    return(list(IAB.obs=IAB.,IAB.exp=mean(MAB), P.attract=P.attract,P.avoid=P.avoid))
    
  } else {
    return(IAB.df)
  }
}



