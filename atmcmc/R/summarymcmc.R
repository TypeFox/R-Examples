#' summary of output from `atmcmc'
#' 
#' Shows a summary of output from an `atmcmc' run. It includes project name, sample mean, functional sample mean, acceptance rate, runtime, and iteration numbers of end of all phases.
#' @param output output from an `atmcmc' run
#' @param name name of the project
#' @return Displays the following items:
#'  \item{name}{name of the project}
#'  \item{estimates}{sample mean}
#'  \item{estimates_of_functional}{functional sample mean}
#'  \item{acceptance_rate}{acceptance rate of the 2nd half of sampling phase}
#'  \item{time_elapsed}{runtime of the full MCMC}
#'  \item{phase_length}{string of iteration numbers to show when each phase has ended. For the sampling phase, this shows when the 1st half of sampling phase has ended also}
#' @examples ## see examples in `atmcmc'
#' @export
summarymcmc<-function(output, name='MCMC'){
  summary = list(name=name, estimates = output$means, estimates_of_functional = output$functionalmeans, acceptance_rate = output$acceptrate,
                 time_elapsed=output$runtime, phase_length=output$sumchain)
  return(summary)
}