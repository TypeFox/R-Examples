liquidlimit <-
function(theta, n)
{
   if (!inherits(theta, "numeric") || theta < 0 || theta > 1)
      stop ("'theta' must be a value between 0 and 1!")
   out <- theta * (n / 25) ^ 0.12
   return(out)
}
