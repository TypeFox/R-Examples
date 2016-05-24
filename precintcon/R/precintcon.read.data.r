#' @export
precintcon.read.data <- function(file, sep = ",", dec = ".", header = TRUE, na.value = NA) {

	data <- read.table(file, header=header, sep=sep, dec=dec)
	
	if (!(length(na.value) == 1 && is.na(na.value)))
	  for (i in 1:length(na.value))
	    for (j in 1:ncol(data))
	      data[which(data[j] == na.value[i]),j] <- NA
	
	data <- data[!is.na(data[1]),]
	data <- data[!is.na(data[2]),]
	
	v <- table(data[1])

	if (length(v[v != 12]) > 0) {
		
		v <- row.names(as.data.frame(v[v != 12]))

		m <- ""

		if (length(v) > 1)
			m <- "Inconsistent data for the years:"
		else
			m <- "Inconsistent data for the year:"

		for (i in 1:length(v))
			m <- paste(m, v[i])
		
		stop(m)
	}
	
	cnames <- c("year", "month")
	
	if (ncol(data) >= 33) {
		
	  data <- data[,1:33]
	  
	  for(i in 1:31)
	    cnames <- c(cnames, paste("d", i, sep = ""))
	  
		class(data) <- c(class(data), "precintcon.daily")		
		
	} else if (ncol(data) >= 3) {
		
	  data <- data[,1:3]
	  
		class(data) <- c(class(data), "precintcon.monthly")	
		
	} else {
		
		stop("Invalid data. Please, check your input file. It should has 3 or 33 columns.")
	}
	
	colnames(data) <- cnames
	
	for (i in 1:ncol(data))
	  data[grep(" +", data[,i]),i] <- NA
	
	return(data)
}
