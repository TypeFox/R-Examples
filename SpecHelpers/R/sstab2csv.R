sstab2csv <-function(in.file = "", out.file = "") {
	
	df <- read.table(in.file, skip = 17, sep = "\t", nrows = 2048)
	write.table(df, file = out.file, row.names = FALSE,
		col.names = FALSE, quote = FALSE, sep = ",")

	}

