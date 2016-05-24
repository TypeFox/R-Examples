#' Set the default options for the MCMC simulation.
#' 
#' This function creates a list with the default options of the MCMC
#' simulation.
#' 
#' 
#' @return A list of default options with elements: \item{lmax}{Maximum number
#' of parent nodes. Default=5.} \item{m}{Number of repeated measurements.
#' Default=1 (no repeats).} \item{dyn}{Lag for the DBN model. Default = 1 when
#' X(t) depends on the previous measurement X(t-1), but dyn can be chosen equal
#' to 2, 3, ...} \item{minPhase}{Minimal length of a segment. Default=2.}
#' \item{maxCP}{Maximal number of changepoints. Default=10.}
#' \item{maxTF}{Maximal number of incoming edges for each node. Default=5.}
#' \item{alphaCP}{Hyperparameter for the number of changepoints. Default=1.}
#' \item{betaCP}{Hyperparameter for the number of changepoints. Default=0.5.}
#' \item{alphaTF}{Hyperparameter for the number of incoming edges. Default=1.}
#' \item{betaTF}{Hyperparameter for the number of incoming edges. Default=0.5.}
#' \item{burnin}{Whether to include a burnin period. Default=F.}
#' \item{psrf.check}{Whether to calculate the potential scale reduction factor
#' (PSRF). Default=F.} \item{pp.l1}{Proposal frequency for level-1
#' hyperparameter moves. Default=0.2.} \item{pp.l2}{Proposal frequency for
#' level-2 hyperparameter moves. Default=0.01.} \item{save.by.node}{Whether to
#' save results separately for each target node. Default=F.}
#' \item{save.file}{Whether to save the results to a file. Default=F.}
#' \item{hyper.fixed}{Whether to keep the network structure prior
#' hyperparameters fixed. Default=F.} \item{cp.fixed}{Whether to keep the
#' changepoints fixed. Default=F.} \item{hyper.init}{Initial values for the
#' network structure prior hyperparameters. Default=NULL.}
#' \item{cp.init}{Initial values for the changepoint locations. Default=NULL.}
#' @author Frank Dondelinger
#' @examples
#' 
#' # Set options to allow saving network and changepoint samples to file
#' options = defaultOptions()
#' options$save.file = TRUE
#' 
#' # NOT EXECUTED
#' # result.bino2 = EDISON.run(dataset$sim_data, 
#' #                  information.sharing='bino_hard',
#' #                  num.iter=5000, output.file='bino2.results',
#' #                  options=options)
#' 
#' @export defaultOptions
defaultOptions <-
function() {
  
  # Maximum number of parent nodes (fan-in restriction). 
  lmax = 5;
  
  # Number of repeated measurements (=1 when no repetition)
  m = 1
  
  # Lag for the DBN model (usually = 1 when X(t) depends on the previous measurement X(t-1), but dyn can be chosen equal to 2, 3, ... )
  dyn = 1 
  
  # Minimal length of a segment (or a phase)
  minPhase = 2
  
  # Maximal number of CPs
  maxCP = 10
  
  # Maximal number of incoming edges (TF) for each node 
  maxTF = 5
  
  # Hyperparameters for the number of CP and incoming edges (TF) 
  alphaCP = 1
  betaCP = 0.5
  alphaTF = 1
  betaTF = 0.5
  
  # By default, do not keep betas constant during burnin
  burnin = FALSE
  
  # Whether to calculate the potential scale reduction factor (PSRF)
  psrf.check = FALSE
  
  # Proposal probabilities for hyperparameter moves
  pp.l2 = 0.01
  pp.l1 = 0.2

  # Output options
  save.by.node = FALSE
  save.file = FALSE

  hyper.fixed = FALSE
  cp.fixed = FALSE

  hyper.init = NULL
  cp.init = NULL

  return(list(lmax=lmax, m=m, dyn=dyn, minPhase=minPhase,
              maxCP=maxCP, maxTF=maxTF, alphaCP=alphaCP, betaCP=betaCP,
              alphaTF=alphaTF, betaTF=betaTF, burnin=burnin,
              psrf.check=psrf.check, save.by.node=save.by.node, save.file=save.file,
              hyper.fixed=hyper.fixed, hyper.init=hyper.init, cp.fixed=cp.fixed,
              cp.init=cp.init, pp.l1=pp.l1, pp.l2=pp.l2))
}

