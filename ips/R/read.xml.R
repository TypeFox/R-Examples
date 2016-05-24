read.xml <- function(x){
	x <- scan(x, what = "c", sep = "\n")
	start <- grep("<alignment id=", x)
	stop <- grep("</alignment>", x)
	
	x <- x[(start + 2):(stop - 2)]
	x <- gsub("\t", "", x)
	x <- x[-grep("sequence", x)]
	
	ind <- seq(1, length(x) - 1, by = 2)
	nb <- length(ind)
	nam <- x[ind]
	nam <- gsub("<taxon idref=\"|\"/>", "", nam)
	seq <- x[ind + 1]
	x <- list(nb = nb, seq = seq, nam = nam)
	class(x)  <- "alignment"
	x <- as.DNAbin(x)
	x
}

