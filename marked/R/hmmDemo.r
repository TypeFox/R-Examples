
#' HMM computation demo functions
#' 
#' Uses fitted hmm model to construct HMM state vectors alpha and phi for demonstration purposes
#' @param object fitted hmm model
#' @param ddl design dat list; if NULL it is created
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @param obs.names names for observations used to label output; if NULL uses ObsLevels
#' @return hmm demo list which includes 1) lnl - the log-likelihood value, 2) alpha - forward probabilities,
#' 3) beta - backward probabilities, 4) phi - scaled forward probabilities, 5) v- intermediate calculation for phi,
#' 6) dmat - 3-d array with observation probability matrix for each occasion, 7) gamma - 3-d array with state transition probability
#' matrix for each occasion, 8) stateprob - predicted state probabilities, 9) local_decode - state predictions for each occasion and individual, 
#' 10) global_decode - state predictions for entire sequence for each individual.
#' @author Jeff Laake
#' @export hmmDemo
#' @useDynLib marked
#' @keywords models
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # cormack-jolly-seber model
#' data(dipper)
#' mod=crm(dipper,model="hmmcjs")
#' x=hmmDemo(mod,state.names=c("Alive","Dead"),obs.names=c("Missed","Seen"))
#' par(mfrow=c(2,1))
#' barplot(t(x$alpha[45,,]),beside=TRUE,names.arg=x$chforwardstrings)
#' barplot(t(x$phi[45,,]),beside=TRUE,names.arg=x$chforwardstrings)
#' # multi-state example showing state predictions
#' data(mstrata)
#' mod=crm(mstrata,model="hmmMSCJS",strata.labels=c("A","B","C"))
#' #' x=hmmDemo(mod)
#' # state predictions are normalized by likelihood value which = rowSums(alpha*beta)
#' cat(paste("\nrowsums = ",rowSums(x$alpha[45,,]*x$beta[45,,],na.rm=TRUE)[2],
#'    "which matches likelihood value",exp(x$lnl[45]),"\n"))
#' # state predictions given the data
#' x$stateprob[45,,]
#' }
hmmDemo <- function(object,ddl=NULL,state.names=NULL,obs.names=NULL)
{
	result=loglikelihood(object,ddl=ddl)
	result$beta=backward_prob(object,ddl=ddl)
	result$stateprob=result$alpha*result$beta/exp(result$lnl)
	result$local_decode=local_decode(object,ddl=ddl)
	result$global_decode=global_decode(object,ddl=ddl)
	if(is.null(state.names))
		if(!is.null(object$data$strata.labels))
			state.names=c(object$data$strata.labels,"Dead")
		else
			state.names=c("Alive","Dead")
	if(is.null(obs.names))obs.names=object$data$ObsLevels
    dimnames(result$alpha)[3]=list(state.names)
	dimnames(result$phi)=dimnames(result$alpha)
	dimnames(result$v)=dimnames(result$alpha)
	dimnames(result$beta)=dimnames(result$alpha)
	dimnames(result$stateprob)=dimnames(result$alpha)
	dimnames(result$gamma)[3:4]=c(list(state.names),list(state.names))
	names(dimnames(result$gamma))=c("Id","Occasion","From_state","To_state")
	dimnames(result$dmat)[3:4]=c(list(obs.names),list(state.names))
	names(dimnames(result$dmat))=c("Id","Occasion","Observation","State")
	names(dimnames(result$alpha))=c("Id","Occasion","State")
	names(dimnames(result$phi))=names(dimnames(result$alpha))
	names(dimnames(result$beta))=names(dimnames(result$alpha))
	names(dimnames(result$stateprob))=names(dimnames(result$alpha))
	names(dimnames(result$v))=names(dimnames(result$alpha))
#	ch=splitCH(object$data$data$ch)
#	result$chfowardstrings=apply(ch,1,function(z) sapply(1:length(z),function(x)paste(z[1:x],collapse="")))
#	result$chbackwardstrings==apply(ch,1,function(z) sapply(1:length(z),function(x)paste(rev(z)[1:x],collapse="")))
	return(result)
}
