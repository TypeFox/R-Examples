
Rhpc_serialize <- function(obj)
{
  .Call("Rhpc_serialize", obj, PACKAGE="Rhpc")
}

Rhpc_unserialize <- function(obj)
{
  .Call("Rhpc_unserialize", obj, PACKAGE="Rhpc")
}

