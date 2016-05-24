update.CADFtest <- function(object, change, ...)
{
# object: is an object of class CADFtest
# change: list of character. It is the change to the model. Ex: list("+x3 -x2", "kernel='Parzen')
	nc <- nchar(object$data.name)
    if (substring(object$data.name, nc, nc)==".") substring(object$data.name, nc, nc) <- "1"
    for (i in 1:length(change))
	{
		if ((substring(change[i], 1, 1)=="+")|(substring(change[i], 1, 1)=="-"))
		{
			# change formula
			newformula <- update.formula(object$data.name, paste("~ .",change[i]))
		    object$call$model <- newformula
		}
		else
		{
			text <- paste("object$call$", change[i], sep="")
			eval(parse(text=text))
		}
	}
	eval(object$call)
}

	
