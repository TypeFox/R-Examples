bigRR <-
function(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, weight = NULL,
		family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
		tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) UseMethod("bigRR")

.onAttach <- 
		function(lib, pkg, ...)
{
	pkgDescription <- packageDescription(pkg)
	pkgVersion <- pkgDescription$Version
	pkgDate <- pkgDescription$Date
	pkgName <- pkgDescription$Package
	pkgTitle <- pkgDescription$Title
	pkgAuthor <- pkgDescription$Author
	pkgMaintainer <- pkgDescription$Maintainer
	packageStartupMessage(paste("\n", pkgName, ": ", pkgTitle, sep = ""))
	packageStartupMessage(paste("Version ", pkgVersion, " (", pkgDate, ") installed", sep = ""))
	packageStartupMessage(paste("Authors: ", pkgAuthor, sep = ""))
	packageStartupMessage(paste("Maintainer: ", pkgMaintainer, "\n", sep = ""))
	packageStartupMessage('Use citation("bigRR") to know how to cite our work.\n')
	packageStartupMessage('!! NOTE !! The bigRR.update() function in bigRR <= 1.3-4 is now bigRR_update().')
	packageStartupMessage('           Please replace in all your old source code.')
	packageStartupMessage('!! NOTE !! The GPU option is only maintained in the R-Forge versions of "bigRR".\n')
}