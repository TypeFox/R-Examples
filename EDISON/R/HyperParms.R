#' Sets up initial values of hyperparameters.
#' 
#' This function initialises the variable HYPERvar with values for the various
#' hyperparameters in the model.
#' 
#' 
#' @param options MCMC settings, possibly from \code{\link{defaultOptions}}.
#' @return Settings for the HYPERvar variable: \item{cD}{Proportion of
#' changepoint moves proposed.} \item{alphaD}{Prior settings for the number of
#' changepoints.} \item{betaD}{Prior settings for the number of changepoints.}
#' \item{c}{Ratio of changepoint birth/death moves proposed.} \item{v0}{Prior
#' settings for the sigma squared parameters.} \item{gamma0}{Prior settings for
#' the sigma squared parameters.} \item{alphad2}{Prior settings for the
#' signal-to-noise ratio delta squared.} \item{betad2}{Prior settings for the
#' signal-to-noise ratio delta squared.} \item{alphalbd}{Prior settings for the
#' number of transcription factors.} \item{betalbd}{Prior settings for the
#' number of transcription factors.}
#' @author Sophie Lebre
#' 
#' Frank Dondelinger
#' @export HyperParms
HyperParms <-
function(options){

  alphaCP = options$alphaCP; betaCP = options$betaCP
  alphaTF = options$alphaTF; betaTF = options$betaTF
  dyn = options$dyn
  smax = options$maxCP; qmax = options$maxTF
  
  #####################################################################
  ### rjMCMC hyperparameters 
  #####################################################################

  ### Level 1 (4 moves: CP birth, CP death, CP update or phase update)

  ## For birth/death/move/update phase acceptance rates
  
  # Unless changepoint sampling is enabled, set proportion of 
  # CP moves proposed to 0
  if(options$cp.fixed) {
    cD = 0 
  } else {
    cD = 0.1
  }
  
  ### level 2 (4 moves)(for each current phase: Pred birth, Pred death or Regression Coefficient update)
  ## for Pred birth/death acceptation
  c = 0.5

  ## For each hidden state (model selection)
  # sig2 ~ IG (v0/2,gamma0/2)
  v0 = 1
  gamma0 = 0.1
  
  # For the signal-to-noise ratio
  # delta2 ~ IG(alphad2,betad2)
  alphad2 = 2
  betad2 = 5
  #######################################################################
  #######################################################################
  
 
  ## For the number of changepoints (CP)
  # for D sampling (D ~ Ga(alphaD,betaD)
  alphaD = alphaCP
  betaD = betaCP
  
  ## For the number of Transcription Factors (TF)
  ## for lambda ~ Ga(alphalbd,betalbd)
  alphalbd = alphaTF
  betalbd = betaTF

  HYPERvar = list(cD=cD, alphaD=alphaD, betaD=betaD, c=c, v0=v0, 
                  gamma0=gamma0, alphad2=alphad2, betad2=betad2, 
                  alphalbd=alphalbd, betalbd=betalbd)
  
  return(HYPERvar)
}

