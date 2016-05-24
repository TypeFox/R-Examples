#' sighsmm_model class
#'
#' A stochastic model of signal peptide produced by \code{signalHsmm}.
#'
#' @details Always a named list of five elements
#' \enumerate{
#' \item \code{aa_group} encoding of amino acids. See \code{\link{aaaggregation}} 
#' for an example.
#' \item \code{pipar} probabilities of initial state in Markov Model.
#' \item \code{tpmpar} matrix of transition probabilities between states.
#' \item \code{od} matrix of response probabilities. Eg. od[1,2] is a probability 
#' of signal 2 in state 1.
#' \item \code{overall_probs_log} probabilities of amino acids in mature protein.
#' \item \code{params} matrix of probability distribution for duration. Eg. params[10,2] 
#' is probability of duration of time 10 in state 2.
#' }
#' @seealso \code{\link{train_hsmm}} \code{\link{predict.sighsmm_model}}
#' @name hsmm_pred
NULL