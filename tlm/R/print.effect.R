print.effect <-
function(x, ...)
 {
  if (!inherits(x, "effect"))
     stop("argument 'x' must be of class 'effecttlm'")
  cat("\n")
  if (attr(x, "unusegarg") == TRUE)
     cat("WARNING: any effect corresponds to the provided arguments.\n\n")
  cat(x$info)
  print(x$effect, ...)
  if (attr(x, "modeltype") > 0)
   {
   	cat("\n")
   	cat("For further information on interpreting the effect use effectInfo().\n\n")
   }
 }
