"view" <-

function(x, maxlines = 10.)



{



	if(NROW(x) <= maxlines) x



	else if(is.matrix(x) || is.data.frame(x))

	

	{



		cat(NROW(x) - maxlines,"remaining rows have not been listed", fill = TRUE)



		x[1:maxlines, ]	



	}



	else if(is.list(x) || is.vector(x))

	

	{



		cat(NROW(x) - maxlines,"remaining rows have not been listed", fill = TRUE)



		x[1:maxlines]	



	}



	else x



}

