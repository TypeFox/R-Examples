
###################################################################################################
#                                                                                                 #
#                                         print.qaw                                               #
#                                                                                                 #
###################################################################################################

print.qaw <- function(x, ...)
{
  cat("\n")
  cat("\n$socket: ")
  print(x$socket)
  cat("\n$banque: ")
  #cat(get("bankName",env=.GlobalEnv)) # Ca pas bon
  cat("\n$call: ")
  print(x$call)
  cat("$name: ")
  print(x$name)
  cat("\n")
  sumry <- array("", c(1, 4), list(1, c("list", "length", "mode", "content")))
  sumry[1, ] <- c("$req",length(x$req),"character","sequences")
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
}


