"listString" <-
function(x, period = FALSE, verbose = FALSE)
{
   if(verbose) cat("\n      entering listString\n")
   flush.console()

   if(!is.character(x))  x <- as.character(x)

   numElements <- length(x)
   out <- if(length(x) > 0)
   {
      switch(
         min(numElements, 3),
         x,
         paste(x, collapse = " and "),
         {
            x <- paste(
               x,
               c(
                  rep(",", numElements - 2),
                  " and",
                  ""),
               sep = "")
            paste(x, collapse = " ")
         })
   } else ""

   if(period) out <- paste(out, ".", sep = "")
   if(verbose) cat("      leaving  listString\n\n")
   flush.console()
   out
}
