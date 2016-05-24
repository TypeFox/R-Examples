#' Bootstrapped Recommendation
#' 
#' Provides recommendations for consistency score and configurational n thresholds to attain a desired level of confidence in a QCA algorithm application.
#' @import QCAGUI bootstrap
#' @importFrom graphics hist
#' @importFrom stats glm plogis predict quantile
#' @importFrom utils flush.console
#' @param qca.data the QCA data frame.
#' @param outcome the outcome variable in the QCA data frame of causal conditions; \code{"OUT"} is the outcome variable for an application of QCA.
#' @param type of QCA application, \code{"crisp"} or \code{"fuzzy"} sets. Default set to \code{type = "crisp"}.
#' @param inclcut range of consistency scores for inclusion. If not specified, this defaults to \code{seq(from = 0.5, to = 1, by = 0.01)}.
#' @param ncut configurational n levels to simulate. Can be altered to give options for the range of minimum to maximum \code{ncut} value that the truth table yields, by naming the the truth table object (e.g. \code{truth}) and calling the minimum and maximum number of cases, using \code{ncut=min(truth$tt$n):max(truth$tt$n)} identified by the truth table. Default set to \code{ncut=2}.
#' @param neg.out [from QCAGUI package] ``Logical, use negation of outcome (ignored if data is a truth table object).'' Default set to \code{neg.out=F}.
#' @param sim number of simulations to run for each combination of parameters. The final number of simulations is \code{length(inclcut)*length(ncut)*sim*2}. Default set to \code{sim=10}.
#' @param verbose prints the system time used to run the simulation and the percent complete. Default set to \code{verbose=T}.
#' @return Significance levels reached (.10,.05, .01, .001) when specifying a combination of inclcut, ncut, and neg.out in a QCA model.
#' @examples
#' data(rallies)
#' P<-rallies$P
#' R<-rallies$R
#' C<-rallies$C
#' U<-rallies$U
#' 
#' qca.data<-data.frame(P,R,C,U)
#' truth<-truthTable(qca.data,outcome="P",sort.by="incl",incl.cut1=0.7,show.cases=TRUE)
#' truth
#' mod1 <- eqmcc(truth,details=TRUE,show.cases=TRUE)
#' mod1
#' 
#' brQCA(qca.data,outcome="P",ncut=1,sim=1)
#' @export
brQCA<-function(qca.data, outcome="OUT", type="crisp", inclcut = "", ncut=2, neg.out=F, sim=10, verbose=T){
  s.data<-sim.brQCA(qca.data, outcome, inclcut = inclcut, ncut=ncut, sim=sim, neg.out=F, type=type, verbose=verbose)
  results<-conf.table(s.data, ncut)
  return(results)
}