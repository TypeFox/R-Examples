#<<BEGIN>>
evalmcmod <- function(expr, nsv = ndvar(), nsu = ndunc(), seed = NULL)
#TITLE Evaluates a Monte-Carlo model
#DESCRIPTION
# Evaluates a \code{\link{mcmodel}} object (or a valid expression) using a specified number of simulations and with (or without) a specified seed.
#KEYWORDS methods
#INPUTS
#{expr}<<A model of class \code{\link{mcmodel}} or a valid expression.>>
#[INPUTS]
#{nsv}<<The number of simulations in the dimension of variability used in the evaluation.>>
#{nsu}<<The number of simulations in the dimension of uncertainty used in the evaluation.>>
#{seed}<<The random seed used for the evaluation. If \samp{NULL} the \samp{seed} is unchanged.>>
#VALUE
# The results of the evaluation. It should be a \samp{mc} object.
#DETAILS
# The model is evaluated. The intermediate variables used to build the \samp{mc} object are not stored.</>
#NOTE
#The seed is set at the beginning of the evaluation. Thus, the complete similarity
#of two evaluations with similar seed is not certain, depending on the structure of your model.
#SEE ALSO
#\code{\link{mcmodel}}</>
#\code{\link{evalmccut}} to evaluate high dimension Monte Carlo Model in a loop.
#EXAMPLE
#data(ec)
#ec$modEC1
#evalmcmod(ec$modEC1,nsv=100,nsu=100,seed=666)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
  if(!is.null(seed)) set.seed(seed)
  Oldv <- ndvar()
  Oldu <- ndunc()
  ndvar(nsv)
  ndunc(nsu)
  x <- try({
	  x <- eval(expr)
	  if(!is.mc(x)) stop("expr does not lead to a mc")
    x}, silent=TRUE)

  ndvar(Oldv)
  ndunc(Oldu)
  if(inherits(x,"try-error")) stop(x,call. = FALSE)
  return(x)
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

