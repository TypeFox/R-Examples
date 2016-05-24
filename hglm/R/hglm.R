`hglm` <-
		function(X = NULL, y = NULL, Z = NULL, family = gaussian(link = identity),
				rand.family = gaussian(link = identity), method = "EQL", conv = 1e-6, maxit = 50, 
				startval = NULL, fixed = NULL, random = NULL, X.disp = NULL, disp = NULL, 
				link.disp = "log", X.rand.disp = NULL, rand.disp = NULL, link.rand.disp = "log", 
				data = NULL, weights = NULL, fix.disp = NULL, offset = NULL, RandC = ncol(Z), 
				sparse = TRUE, vcovmat = FALSE, calc.like = FALSE, bigRR = FALSE, verbose = FALSE, ...) UseMethod("hglm")

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
	#CRANpkg <- try(available.packages(contriburl = 'http://cran.at.r-project.org/bin/macosx/contrib/3.0'), silent = TRUE)
	#if (!is.null(CRANpkg) & class(CRANpkg) != "try-error" & nrow(CRANpkg) != 0) {
	#	cranVersion <- CRANpkg[pkgName, 'Version']
	#	if (pkgVersion != cranVersion) {
	#		packageStartupMessage(paste(
	#					"The installed ", pkgName," version (", pkgVersion, ") is not the same as the stable\n",
	#					"version available from CRAN (", cranVersion, "). Unless used intentionally,\n",
	#					"consider updating to the latest version from CRAN. For that, use\n",
	#					"'install.packages(\"", pkgName, "\")', or ask your system administrator\n",
	#					"to update the package.\n", sep = ""))
	#	}
	#}
	packageStartupMessage('Use citation("hglm") to know how to cite our work.\n')
	packageStartupMessage('Discussion: https://r-forge.r-project.org/forum/?group_id=558')
	packageStartupMessage('BugReports: https://r-forge.r-project.org/tracker/?group_id=558')
	packageStartupMessage('VideoTutorials: http://www.youtube.com/playlist?list=PLn1OmZECD-n15vnYzvJDy5GxjNpVV5Jr8')
}

