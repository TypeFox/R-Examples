".onAttach" <- 
function(libname, pkgname) {
	# Attach JDBC Driver JARs
	jarfiles <- list.files(system.file("drivers", package=pkgname))
	jarfiles <- setdiff(jarfiles, c("README", "README.TXT")) # omit READMEs
	if (length(jarfiles) > 0) {
		packageStartupMessage("Attaching jar files in \"drivers\" folder: \n", paste("\t", jarfiles, "\n"))
		loadJDBCDriver(system.file("drivers", jarfiles, package=pkgname))
	}
	invisible()	
}

