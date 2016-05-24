.centr <- function(X,m) t (t(X) - m)


.mymean <- function (X)
{
	if (any (!is.na (X)))
		mean (X[!is.na (X)], trim = 0.4)
	else
		1
}

#.scale <- function (X)
#{
#	med <- apply (X, 2, median)
#	dev <- abs (.centr (X, med))
#	dev[dev == 0] <- NA	
#	apply(dev, 2, .mymean)
#}

# .First.lib <-
# function(lib,pkg)
# {	##	2do: delete, using namespaces now; this function is useless.
#    library.dynam("pcaPP",pkg,lib)
#    library(mvtnorm)
#    cat("pcaPP 0.1-1 loaded\n")
# }

.ParseControlStructure = function (control, arguments)
{
	if (!is.list(control))
		stop ("Invalid argument type: control structure must be of type list")

	if (missing(arguments))
		arguments = attributes(control)$names

	for (curname in arguments)
	{
		if (!is.null (control[[curname]]))	##	 if this argument is provided in the control structure:
			eval (parse (text = paste ("eval.parent (substitute(", curname, "<- control$", curname, "), n = 3)", sep = "")))
	}
}

.ParseDevString = function (method)
{
	if (method[1] == "mad")  return (0)
	if (method[1] == "sd")  return (1)
	if (method[1] == "Qn" | method[1] == "qn" )  return (2)
	return (1)
}

.Conv2Matrix <- function (x, sx = substitute (x))
{
	if(is.matrix(x))
		return (x)
	if(is.data.frame(x))
		return (data.matrix(x))
	return (matrix(x, nrow = length(x), ncol = 1, dimnames = list(names(x), deparse(sx))))
}


.colMedians <- function (x)
{
	if (is.null (dim (x)))
		return (median (x))
	apply (x, 2, median)
}

.GetFunctionName <- function (f, ...)
{
	form <- formals (f)
	if (!is.null (form$fa.Name))
		return (eval (form$fa.Name))

	if (!is.null (form$NAME))
		return (f (NAME = TRUE, ...))

	if (!is.null (attributes (f)$NAME))
		return (attributes (f)$NAME)

	return (NULL)
}

.flush.cat <- function (...)
{
	cat (...)
	flush.console ()	
}



.colSds <- function (x)	##	function for replacing sd until the R-people remove the sd matrix warning message
{
    if (is.data.frame(x)) 
        x <- as.matrix(x)

    if (!is.array (x) || length (dim (x)) != 2)
	stop ("'x' must be an array of two dimensions")

    apply (x, 2, sd)
}
