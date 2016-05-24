#' @export
precintcon.spi.per.year.analysis <- function(
   object, 
   period       = 3, 
   distribution = "Gamma", 
   FUN          = mean
) {

   m <- match.call()
   
   fun.name <- ifelse(is.null(m$FUN), "mean", as.character(m$FUN))
   
   result <- precintcon.spi.analysis(object = object, period = period, distribution = distribution)
   
   result <- aggregate(result$spi, by = list(result$year), FUN = FUN)

   colnames(result) <- c("year", paste("spi", fun.name, sep = "."))
   
   class(result) <- c("data.frame", "precintcon.spi.per.year")

   return(result)
}
