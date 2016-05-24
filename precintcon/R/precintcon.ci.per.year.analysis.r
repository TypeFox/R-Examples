#' @export 
precintcon.ci.per.year.analysis <- function(object, interval = 1.0) {

	if (is.element("precintcon.daily", class(object))) {

		d <- object;
		
		if (nrow(d) %% 12 == 0) {
			
			y <- c();
	
			for (j in seq(1, nrow(d), by=12)) {

				data.fd      <- precintcon.fd(precintcon.classification(d[j:(j+11),], interval));
				data.ci      <- precintcon.ci.analysis(data.fd);
				
				y <- c(y, data.ci$ci);
			}
			
			tmp <- data.frame(unique(d[,1]), y)
			
			colnames(tmp) <- c("year", "ci")

			data <- tmp
			
		} else 
			stop("Invalid data serie. Please, check the number of months in your data. It should be multiple of 12.")
   }
	
	return(data);
}