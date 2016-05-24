#<<BEGIN>>
mcmodel <- function(x, is.expr=FALSE)
#TITLE Monte Carlo model
#DESCRIPTION
# Specify a \samp{mcmodel}, without evaluating it, for a further evaluation using \code{\link{evalmcmod}}.
#KEYWORDS methods
#INPUTS
#{x}<<An \R call or an expression.>>
#[INPUTS]
#{is.expr}<< \samp{FALSE} to send a call,  \samp{TRUE} to send an expression (see Examples)>>
#VALUE
# an \R expression, with class \samp{mcmodel}
#DETAILS
# The model should be put between \samp{\{} and the last line should be of the form \samp{mc(...)}.
# Any reference to the number of simulation in the
#dimension of variability should be done via \samp{ndvar()} or (preferred) \samp{nsv}.
# Any reference to the number of simulations in the dimension of uncertainty
# should be done via \samp{ndunc()} or (preferred) \samp{nsu}.
#SEE ALSO
#\code{\link{expression}}.</>
#\code{\link{evalmcmod}} to evaluate the model.</>
#\code{\link{mcmodelcut}} to evaluate high Dimension Monte Carlo Model in a loop.
#EXAMPLE
#modEC1 <- mcmodel({
#    conc <- mcdata(10,"0")
#    cook <- mcstoc(rempiricalD,values=c(0,1/5,1/50),prob=c(0.027,0.373,0.600))
#    serving <- mcstoc(rgamma,shape=3.93,rate=0.0806)
#    expo <- conc * cook * serving
#    dose <- mcstoc(rpois, lambda=expo)
#    risk <- 1-(1-0.001)^dose
#   mc(conc,cook,serving,expo,dose,risk)
#   })
#evalmcmod(modEC1,nsv=100,nsu=100)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{

  if(!is.expr) x <- as.expression(substitute(x))
  if(!is.expression(x)) stop("x can not be evaluate as an expression")

    last <- x[[1]][length(x[[1]])]
    lastcall1 <- substr(deparse(last,width.cutoff = 500), 1, 3)
    if (lastcall1 != "mc(") warning("The last call should be 'mc(...)'")
    class(x) <- "mcmodel"
    return(x)
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

