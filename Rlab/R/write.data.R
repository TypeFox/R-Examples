"write.data" <-

function(x, file)



{



	length.old <- options()$length



	options(length = 10000.)



	if(missing(file)) {



		file <- paste(as.character(substitute(x)), ".output",



			sep = "")



		cat(" The data set is being written to the UNIX file: ",



			file, fill = TRUE)



	}



	sink(file)



	print(x, quote = FALSE)



	sink()



	options(length = length.old)



	invisible()



}

