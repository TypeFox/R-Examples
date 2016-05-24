#' Set the ordinary differential equation (ODE) for LM-OSL measurements
#' in the enery-band-model of quartz.
#'
#' This function provides the differential equations for LM-OSL measurements. This function
#' is necessary because the stimulation power (P) is time-dependent.
#'
#' With RLumModel, the number of coupled
#' differential equations will automatically be adjusted, because of an identifying
#' feature for electron traps and hole centres: This function identifies if B (see
#' \code{\link{set_Pars}}) is 0 (electron trap) or not (hole trap). For advanced users it is thus
#' possible to change the available sets or construct own parameter sets with
#' arbitrary numbers of electron traps and hole centres without taking care on
#' coding the right syntax of the ODEs, without changing the complete code
#' and without advanced knowledge of R coding.
#'
#' @param t \code{\link{numeric}} (\bold{required}): timesteps (set from 'simulate_...' functions)
#'
#' @param n \code{\link{numeric}} (\bold{required}): concentration of electron-/holetraps, valence- and conductionband
#' from step before.
#'
#' @param parameters.step \code{\link{list}} (\bold{required}): parameters for every specific
#' calculation has different parameters (heatingrate, pair-production-rate, ...)
#' and this information is given to the ODE.
#'
#' @return This function returns a list with all changes in the concentration of electron-/holetraps, valence and conductionband
#'
#' @note This function calculates the ODE for the energy-band-model for LM-OSL measurements.
#'
#' @section Function version: 0.1.0
#'
#' @author Johannes Friedrich, University of Bayreuth (Germany),
#'
#' @references
#'
#' Bailey, R.M., 2001. Towards a general kinetic model for optically and thermally stimulated
#' luminescence of quartz. Radiation Measurements 33, 17-45.
#'
#' Bailey, R.M., 2002. Simulations of variability in the luminescence characteristics of natural
#' quartz and its implications for estimates of absorbed dose.
#' Radiation Protection Dosimetry 100, 33-38.
#'
#' Bailey, R.M., 2004. Paper I-simulation of dose absorption in quartz over geological timescales
#' and it simplications for the precision and accuracy of optical dating.
#' Radiation Measurements 38, 299-310.
#'
#' Pagonis, V., Chen, R., Wintle, A.G., 2007: Modelling thermal transfer in optically
#' stimulated luminescence of quartz. Journal of Physics D: Applied Physics 40, 998-1006.
#'
#' Pagonis, V., Wintle, A.G., Chen, R., Wang, X.L., 2008. A theoretical model for a new dating protocol
#' for quartz based on thermally transferred OSL (TT-OSL).
#' Radiation Measurements 43, 704-708.
#'
#' @examples
#'
#' #so far no example available
#'
#' @noRd
.set_ODE_LM_OSL <- function(
  t,
  n,
  parameters.step
  ){

  ##============================================================================##
  ## unpack parameters to be used in this function and to keep the ODE code clear
  ##============================================================================##

  N <- parameters.step$parms$N
  E <- parameters.step$parms$E
  s <- parameters.step$parms$s
  A <- parameters.step$parms$A
  B <- parameters.step$parms$B
  Th <- parameters.step$parms$Th
  E_th <- parameters.step$parms$E_th
  k_B <- parameters.step$parms$k_B
  W <- parameters.step$parms$W
  K <- parameters.step$parms$K

  b <- parameters.step$b
  a <- parameters.step$a
  R <- parameters.step$R
  P <- parameters.step$P
  temp <- parameters.step$temp
  ##============================================================================##


  with(as.list(c(n,parameters.step)), {

    dn <- numeric(length(N)+2)

    j <- 0;
    jj <- 0;
    for (i in 1:length(N)){
      if (B[i] == 0)    {      #use recombination propability of recombination centers to identify electron traps, because they had no recombination propability to recomibnation centers from conduction band
        j <- j+1
        jj <- jj+1
        dn[i] <- n[length(N)+1]*(N[i]-n[i])*A[i]-n[i]*P*a*t*Th[i]*exp(-E_th[i]/(k_B*(273+temp+b*t)))-n[i]*s[i]*exp(-E[i]/(k_B*(273+temp+b*t)))
      } else {#calculate recombination centers
        jj <- jj+1
        dn[i] <- n[length(N)+2]*(N[i]-n[i])*A[i]-n[i]*s[i]*exp(-E[i]/(k_B*(273+temp+b*t)))-n[length(N)+1]*n[i]*B[i]
      }

    }

    ## conduction band
    dn[length(N)+1] = R-sum(dn[1:j])-sum(n[length(N)+1]*n[(j+1):jj]*B[(j+1):jj])

    ## make sure if conduction band calculation comes from Bailey or not (see papers for differences between ODEs)
    if (parms$model != "Bailey 2001" && parms$model != "Bailey2004" && parms$model != "Bailey2002")
    {
      dn[length(N)+2] = R-sum(dn[(j+1):jj])-sum(n[length(N)+1]*n[(j+1):jj]*B[(j+1):jj])

    } else {
      dn[length(N)+2] = R-sum(dn[(j+1):jj])          # valence band ODE for Bailey model 2001/2002/2004
    }


    return(list(dn))# return the rate of change

           })


}
