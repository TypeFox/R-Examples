#<<BEGIN>>
ndvar <- function(n)
#NAME mc.control
#TITLE Sets or Gets the Default Number of Simulations.
#DESCRIPTION
# Sets or retrieves the default number of simulations.
#KEYWORDS misc
#INPUTS
#{n}<<Number of simulations.>>
#DETAILS
#\samp{ndvar()} gets and \samp{ndvar(n)} sets the default number of simulation in the 1D simulations
#or the number of simulation in the variability dimension in the 2D simulations.</>
#\samp{ndunc()} gets and \samp{ndunc(n)} sets the number of simulations in the uncertainty dimension
#in the 2D simulations.</>
#\samp{n} is rounded to its ceiling value.</>
#The default values when loaded are 1001 for \samp{ndvar} and 101 for \samp{ndunc}.
#VALUE
#The current value, AFTER modification if \samp{n} is present (!= \samp{options}).
#EXAMPLE
#(oldvar <- ndvar())
#(oldunc <- ndunc())
#mcstoc(runif,type="VU")
#ndvar(12)
#ndunc(21)
#mcstoc(runif,type="VU")
#ndvar(oldvar)
#ndunc(oldunc)

#CREATED 08-01-25
#--------------------------------------------
{
  if(!exists("mc.control", envir=.BaseNamespaceEnv))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.BaseNamespaceEnv )
  x <- get("mc.control", envir=.BaseNamespaceEnv)
   if(!is.list(x) || is.null(x$nsv) || is.null(x$nsu))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.BaseNamespaceEnv )
  if(!missing(n)){
    if (n > 0) x$nsv <- ceiling(n)
        else stop("Invalid n")
    assign("mc.control",x, envir=.BaseNamespaceEnv)}
  return(x$nsv)}

#<<BEGIN>>
ndunc <- function(n)
#ISALIAS ndvar
#--------------------------------------------
{
  if(!exists("mc.control", envir=.BaseNamespaceEnv))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.BaseNamespaceEnv )
  x <- get("mc.control", envir=.BaseNamespaceEnv)
   if(!is.list(x) || is.null(x$nsv) || is.null(x$nsu))
    assign("mc.control",list(nsv=1001,nsu=101),envir=.BaseNamespaceEnv )
  if(!missing(n)){
    if (n > 0) x$nsu <- ceiling(n)
        else stop("Invalid n")
    assign("mc.control",x, envir=.BaseNamespaceEnv)}
  return(x$nsu)}
