odfWeaveControl <- function(
   zipCmd = c("zip -r $$file$$ .", "unzip -o $$file$$"),
   cleanup = ! debug,
   verbose = TRUE,
   debug = FALSE)
{
	# pass ... args to figGen here?
   list(
      zipCmd = zipCmd,
      cleanup = cleanup,
      verbose = verbose,
      debug = debug)
}

