
.getVtext <- function (idx)
{

	stopifnot (length (idx) == 1)

	if (idx == 1)			#	the package name
		return ("pcaPP")
	if (idx == 2)			#	the package version
		return (installed.packages ()[.getVtext (1),"Version"])
	if (idx == 3)			#	the matlab functions
		return (c ("l1median\\_HoCr", "l1median\\_VaZh", "PCAgrid", "PCAproj", "qn", "sPCAgrid"))
	if (idx == 4)			#	the example
	{
		fn <- system.file ("doc", "matlab.example.txt", package = .getVtext (1))
		return (paste (readLines (fn), collapse = "\n"))
	}

	stop ("unkown idx value")
}
