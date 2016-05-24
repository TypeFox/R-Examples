# Clear the workspace of all objects whose names don't start with a full stop
clear <-
function()
{
env <- as.environment(1)
to.go <- ls(env, all.names=FALSE)
continue <- TRUE
while (continue) {
  nxt <- search()[[2]]
  # bit of a botch
  if (substr(nxt, 1, 8)!="package:")
    detach()
  else
    continue <- FALSE
  }
remove(list=to.go, envir=env)
}
