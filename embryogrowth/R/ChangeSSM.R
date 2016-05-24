#' ChangeSSM convert Schoolfield-Sharpe-Magnuson model from 4 to 6 parameters or reverse
#' @title Generate set of parameters for Schoolfield-Sharpe-Magnuson model
#' @author Marc Girondot
#' @return A vector with parameters
#' @param temperatures A vector with incubation temperatures in degrees Celsius
#' @param parameters A vector of parameters for model to be converted (4 or 6 parameters)
#' @param initial.parameters A vector of parameters for initial model model to be fited (4 or 6 parameters)
#' @param ... A control list to be used with optim, see ?optim
#' @description Generate a set of parameters for Schoolfield-Sharpe-Magnuson model
#' @examples
#' \dontrun{
#' data(resultNest_6p)
#' x1 <- resultNest_6p$par
#' data(resultNest_4p)
#' x2 <- resultNest_4p$par
#' temperaturesC <- (200:350)/10
#' s <- ChangeSSM(temperatures=temperaturesC, parameters=x1, initial.parameters=x2)
#' plotR(list(resultNest_6p, resultNest_4p, s), ylim=c(0,0.3), 
#' col=list("black", "red", "green"), lty=list(1,1,1), 
#' legend=list("R function to mimic", "Initial new R function", 
#' "Fitted new R function"), show.box=FALSE)
#' # Other example to fit anchored parameters
#' data(resultNest_4p)
#' x0 <- resultNest_4p$par
#' t <- hist(resultNest_4p, plot=FALSE)
#' x <- c(3.4, 3.6, 5.4, 5.6, 7.6, 7.5, 3.2)
#' names(x) <- seq(from=range(t$temperatures)[1], to=range(t$temperatures)[2], 
#'      length.out=7)
#' newx <- ChangeSSM(temperatures = (200:350)/10, parameters = x0, 
#'        initial.parameters = x, 
#'        control=list(maxit=5000))
#'  # Example on how to generate a set of SSM parameters from anchored parameters
#'  xanchor <- GenerateAnchor(nests=resultNest_4p)
#'  x <- resultNest_4p$par
#'  xanchor["294"] <- 0
#'  xanchor["308"] <- 2.3291035
#'  xprime <- ChangeSSM(parameters = xanchor,
#'                      initial.parameters = x, control=list(maxit=5000))
#'  plotR(result=resultNest_4p, parameters=list(resultNest_4p$par, xprime$par), 
#'        ylim=c(0,0.3), col=c("black", "red"), 
#'        legend=list("Fitted parameters", "Constrainted parameters"))
#' }
#' @export

ChangeSSM <- function(temperatures=(200:350)/10, 
                      parameters=stop("A set of parameters must be supplied"), 
                      initial.parameters=stop("A set of parameters for new model must be supplied"), 
                      ...) {

growth.rate <- .SSM(273.15+temperatures, parameters)[[1]]*1E5

c <- list(...)

if (length(c)==0) {
  s <- optim(par=initial.parameters, fn=.fitSSM, temperatures=temperatures, growth.rate=growth.rate)
} else {
  s <- optim(par=initial.parameters, fn=.fitSSM, temperatures=temperatures, growth.rate=growth.rate, control=c[[1]])
}
return(s)
}

