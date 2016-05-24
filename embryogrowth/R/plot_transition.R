#' plot_transition show fonction used for transition
#' @title Show fonction used for transition
#' @author Marc Girondot
#' @return Nothing
#' @param result A result object
#' @param parameters Set of parameters. If both result and parameters are indicated, parameters have priority.
#' @param sizes The range of possible sizes
#' @param ... Parameters for plot() such as main= or ylim=
#' @description Plot the transition function
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' data(resultNest_4p)
#' # Get a set of parameters without transition
#' x1 <- resultNest_4p$par
#' # Generate a set of parameters with transition
#' x2 <- switch.transition(x1)
#' x2 <- x2[names(x2)!="transition_P"]
#' x2["transition_S"] <- 4
#' pfixed <- c(rK=2.093313, transition_P=20)
#' resultNest_4p_transition <- searchR(parameters=x2, fixed.parameters=pfixed, 
#' temperatures=formated, derivate=dydt.Gompertz, M0=1.7, 
#' test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p_transition)
#' # show the model for smallest size
#' plotR(resultNest_4p_transition, ylim=c(0,0.3))
#' # show the model for larger sizes
#' plotR(resultNest_4p_transition, set.par=2, ylim=c(0,0.3))
#' # plot model for both together
#' plotR(resultNest_4p_transition, set.par=c(1,2), ylim=c(0,0.3), 
#'        col=c("red", "black"), legend=list("Initial", "End"))
#' plot_transition(result=resultNest_4p_transition, las=1, sizes=c(0,40))
#' compare_AIC(one.model=list(resultNest_4p), two.models=list(resultNest_4p_transition))
#' # Note that the model with fitted transition_P is trivial. Embryos grow fast until  
#' # they reach hatchling size and then growth rate becomes null!
#' }
#' @export


plot_transition <-
function(result=NULL, parameters=NULL, sizes=c(0,40), ...) {

	parssm <- c(parameters, result$par, result$fixed.parameters)
	
	if (!is.na(parssm["transition_S"]) & !is.na(parssm["transition_P"])) {
	
		L <- modifyList(list(xlim=sizes), list(...))

		x1 <- sizes[1]
		x2 <- sizes[2]
	
		L <- modifyList(list(x=seq(from=x1, to=x2, by=0.1), y=1/(1+exp((1/parssm["transition_S"]*(parssm["transition_P"]-seq(from=x1, to=x2, by=0.1))))), bty="n", type="l", ylim=c(0,1), xlab="Size", ylab="Proportion of model 2"), L)
	

		a <- do.call(plot, L) 
	
	} else {
		print("Transition parameters are not defined.")
	}

}
