ParseControlStructure = function (control, arguments)
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
