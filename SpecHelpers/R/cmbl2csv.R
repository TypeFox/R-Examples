cmbl2csv <-
function(in.file = "", out.file = "") {
	
	raw.data <- readLines(in.file, warn = FALSE)
	good.data <- NULL
	bad.rows <- NULL

# Mark lines with "Z2" or "<" for removal

		for (i in 1:length(raw.data)) {
			if (identical(grep("Z2|<", raw.data[i]), as.integer(1))) bad.rows <- c(bad.rows, i)
			}
	
# Z1 apparently indicates absorbance less than zero, replace with 0

		for (i in 1:length(raw.data)) {
			if (identical(grep("Z1", raw.data[i]), as.integer(1))) raw.data[i] <- 0.0
			}

	good.data <- as.numeric(raw.data[-bad.rows]) # remove the bad lines
	
# Data is in one long column, split into two matching wavelength with absorbance

	n <- length(good.data)
	df <- cbind(cbind(good.data[1:(n/2)], good.data[((n/2)+1):n]))
	write.table(df, file = out.file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

	}

