#' \code{adapt} is a function that performs adaptations and plans the secondary group sequential trial. The
#' effect size used for planning the secondary trial is a weighted mean between the interim estimate \code{theta} and 
#' the initially assumed estimate \code{delta} (pT\$delta) of the primary trial.
#' 
#' @title Adaptations in group sequential trials
#' @export
#' @param pT object of the \code{class} \code{GSTobj}; primary trial design
#' @param iD interim data; a list with the variables \code{T} and \code{z}; list(T = stage of interim analysis, z = interim z-statistic)
#' @param SF spending function for the secondary trial
#' @param phi parameter of spending function for the secondary trial when SF=3 or 4 (See below)
#' @param cp conditional power
#' @param theta new effect size (default: estimate from interim analysis)
#' @param I2min minimal total information of secondary trial
#' @param I2max maximal total information of secondary trial
#' @param swImax maximal incremental information per stage
#' @param delta initially assumed effect size for the primary trial (default: estimate from primary trial)
#' @param weight weight of \code{theta} when updating the effect size estimate as weighted mean of \code{theta} and \code{delta}
#' @param warn option if warnings should be printed to the screen (default: true)
#' @return \code{adapt} returns an object of the \code{class} \code{GSTobj}; the design of the secondary trial. The adaptation rule is as in the first simulation example of Brannath et al.(2008). If no adaptations are performed, the function returns \code{sT} = NULL. An object of \code{class} \code{GSTobj} is a list containing the following components: \item{sT}{ secondary trial}
#' @references Brannath, W, Mehta, CR, Posch, M (2008) ''Exact confidence bounds following adaptive group sequential tests'', \emph{Biometrics} accepted.
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @details
#' If no adaptation is performed then this indicates that the original plan is kept. In this case \code{sT} is set to \code{NULL}. 
#' If an adaptation is performed \code{sT} is a list which contains the following elements:
#' \tabular{ll}{
#'    \code{K}\tab number of stages\cr
#'    \code{a}\tab lower critical bounds of secondary group sequential design(are currently always set to -8)\cr
#'    \code{b}\tab upper critical bounds of secondary group sequential design\cr
#'    \code{t}\tab vector with cumulative information fractions\cr
#'    \code{al}\tab alpha (type I error rate); equal to the conditional type I error rate of the primary trial\cr
#'    \code{SF}\tab spending function\cr
#'    \code{phi}\tab parameter of spending function when SF=3 or 4 (See below)\cr
#'    \code{alab}\tab alpha-absorbing parameter values of secondary group sequential design\cr
#'    \code{als}\tab alpha-values ''spent'' at each stage of secondary group sequential design\cr
#'    \code{Imax}\tab maximum information number\cr
#'    \code{delta}\tab effect size used for planning the secondary trial\cr
#'    \code{cp}\tab conditional power\cr
#' }
#' A value of \code{SF}=3 is the power family. Here, the spending function is \eqn{t^{\phi}}, 
#' where phi must be greater than 0. A value of \code{SF}=4 is the Hwang-Shih-DeCani family, 
#' with the spending function \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where phi cannot be 0.
#' @seealso \code{\link{GSTobj}}, \code{\link{print.GSTobj}}, \code{\link{plot.GSTobj}}, \code{\link{plan.GST}}
#' @examples
#' ##The following performs an adaptation of the sample size and 
#' ##number of interim analyses after the first stage of the primary trial. 
#'
#' pT=plan.GST(K=3,SF=4,phi=-4,alpha=0.05,delta=6,pow=0.9,compute.alab=TRUE,compute.als=TRUE)
#'
#' iD=list(T=1, z=1.090728)
#'
#' swImax=0.0625
#'
#' I2min=3*swImax
#' I2max=3*swImax
#'
#' sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.8,theta=5,I2min,I2max,swImax)
#'
#'
#' @keywords methods
adapt <- function(pT, iD, SF, phi, cp, theta=iD$z/(pT$t[iD$T]*pT$Imax), I2min, I2max, swImax, delta=pT$delta, weight=1, warn=TRUE)
{    
  cerror <- cer(pT=pT, iD=iD)
  K2min <- ceiling(I2min/swImax)
  K2max <- ceiling(I2max/swImax)
  
  if(!is.null(delta) && (delta != 0) && !is.null(weight)) theta <- (1 - weight)*delta + weight*theta

  if(theta <= 0) {
      cat(paste("No design adaptation performed.\n"))
      NULL
  } else {
      I2 <- (qnorm(cp) - qnorm(cerror))^2 / theta^2
      if(I2 <= swImax) {
          cat("\nI2<=swImax, 1 stage trial.\n")		
          sT <- list(K=1, Imax=I2, SF=SF, phi=phi, t=1, a=-8, b=1-qnorm(cerror), alpha=cerror, compute.alab=FALSE, compute.als=FALSE)
          cpower <- cp
          sT$delta <-theta
      } else{
          K2 <- K2min
          I2 <- K2*swImax
          
          while(K2 <= K2max) {
              sT <- plan.GST(K=K2, SF=SF, phi=phi, alpha=cerror, Imax=I2, compute.alab=FALSE, compute.als=FALSE)
              sT$delta <- theta
              cpower <- cp(sT)
              
              if(cpower > cp) {
                  sT <- plan.GST(K=ifelse(K2>K2min,K2-1,K2), SF=SF, phi=phi, alpha=cerror, delta=theta, pow=cp, compute.alab=FALSE, compute.als=FALSE)	
                  cpower <- cp(sT)
                  break
              } else if(cpower <= cp && K2 >= K2max) {
                  if(warn) warning("Power may be lower than planned.")
                  break
              } else{
                  K2 <- K2+1
                  I2 <- K2*swImax
              }
          }
      }
      sT$cp <- cpower
      class(sT) <- "GSTobj"
      sT
  }
}

