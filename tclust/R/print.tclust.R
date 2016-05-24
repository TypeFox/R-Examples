print.tclust <-
function (x, ...)
{
  cat ("* Results for TCLUST algorithm: *\n")
  cat ("trim = ", x$par$alpha, ", k = ", x$k, "\n", sep = "")

  cat ("Classification (trimmed points are indicated by 0", "):\n")

  print (x$cluster)
  cat ("Means:\n")
  print (x$centers)
  if (x$obj < (-1e+20)) 
    warning ("The solution is not reliable. More iterations are probably needed.")
  cat ("\nTrimmed objective function: ", x$obj, "\n")

  if (!is.null (x$restr.fact))
    cat ("Selected restriction factor:", x$restr.fact, "\n")
  cat (round (x$int$iter.converged / x$int$iter.successful* 100), "% of iterations converged successfully.\n", sep= "")
}
