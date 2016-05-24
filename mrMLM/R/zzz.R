.onLoad <- function(...) {
  ## use .onLoad rather than .First.lib when there is namespace
  #cat("\nUse 'mrMLM()' to restart the programe.\n",fill=TRUE)
  if (interactive()) mrMLM()
}

 


