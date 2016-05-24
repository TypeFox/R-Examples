#' @name bindPosterior
#' @importFrom dplyr bind_rows
#' @export bindPosterior
#' 
#' @title Bind Posterior Distributions
#' @description After determining the posterior distributions are satisfactory,
#'   it can be advantageous to bind the posterior distributions together in
#'   order to aggregate values and perform other manipulations and analyses.
#'   
#' @param hydePost An object of class \code{HydePosterior}
#' @param relabel_factor Logical.  If \code{TRUE}, factors that had been 
#'   converted to integers for the JAGS code can be relabelled as factors 
#'   for additional analysis in R.
#'   
#' @details For the purposes of this function, it is assumed that if the 
#'   posterior distributions are satisfactory, the multiple chains in a run 
#'   can be bound together.  Subsequently, the multiple runs are bound 
#'   together.  Lastly, the factors are relabeled, if requested.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples
#' #' data(PE, package="HydeNet")
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat,
#'                      data = PE) 
#'   
#'                  
#' compiledNet <- compileJagsModel(Net, n.chains=5)
#' 
#' #* Generate the posterior distribution
#' Posterior <- HydePosterior(compiledNet, 
#'                            variable.names = c("d.dimer", "death"), 
#'                            n.iter=1000)
#' 
#' Bound <- bindPosterior(Posterior)
#' 
#' #* Bind a Decision Network
#' #* Note: angio shouldn't really be a decision node.  
#' #*       We use it here for illustration
#' Net <- setDecisionNodes(Net, angio, treat)
#' compiledDecision <- compileDecisionModel(Net, n.chains=5)
#' PosteriorDecision <- HydePosterior(compiledDecision, 
#'                                    variable.names = c("d.dimer", "death"),
#'                                    n.iter = 1000)
#' 
bindPosterior <- function(hydePost, relabel_factor=TRUE){

  if (class(hydePost$codas) == "mcmc.list")
    bound <- dplyr::bind_rows(lapply(seq_along(hydePost$codas), bind_chains_mcmclist, hydePost))
  else 
    bound <- dplyr::bind_rows(lapply(hydePost$codas, bind_chains_list))
  
  if (relabel_factor){
    factors_to_relabel <- names(bound)[names(bound) %in% names(hydePost$factorRef)]
    for(i in factors_to_relabel){
      bound[i] <- factor(bound[[i]], 
                         levels=hydePost$factorRef[[i]]$value,
                         labels=hydePost$factorRef[[i]]$label)
    }
  }

  as.data.frame(bound)
}



#**** UTILITY FUNCTIONS
bind_chains_mcmclist <- function(mcmc, hydePost){
  as.data.frame(hydePost$codas[[mcmc]]) %>%
    dplyr::mutate_(chain_index = ~mcmc,
            obs_index = ~1:n())
}

bind_chains_list <- function(mcmc){
  lapply(1:length(mcmc),
         function(chain)
           as.data.frame(mcmc[[chain]]) %>%
             dplyr::mutate_(chain_index = ~chain,
                     obs_index = ~1:n())) %>%
    dplyr::bind_rows()
}
