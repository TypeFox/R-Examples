Rhpc_enquote<- function(...)
{
  args <- list(...)
  .Call("Rhpc_enquote", args, PACKAGE="Rhpc")
}

Rhpc_splitList <- function(var,num)
{
  .Call("Rhpc_splitList", as.list(var),as.integer(num), PACKAGE="Rhpc")
}

