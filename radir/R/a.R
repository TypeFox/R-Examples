a <-
function(x, pos=1, envir=as.environment(pos))
{
  return(get("u.local2", envir=envir)(x)*b(x))
}
