.onLoad <- function(libname, pkgname) {
	
	.Tcl('set ng_windowManager("ngInstance") 0')
	
	## graph display
	tclfile <- file.path(find.package(package = pkgname, lib.loc = libname),"tcl", "GraphDisplay.tcl")
	tcl("source", tclfile)
	
	## image resizing function in C 
	.Tcl(paste('load "',system.file("libs",.Platform$r_arch,paste("ImgscaleTea",.Platform$dynlib.ext,sep=''), package = pkgname, lib.loc = libname),'"',sep=''))
	.Tcl(paste('load "',system.file("libs",.Platform$r_arch,paste("DisplaystuffTea",.Platform$dynlib.ext,sep=''),package = pkgname, lib.loc = libname),'"',sep=''))
	
	
	## tk2d display
	tclfile <- file.path(find.package(package = pkgname, lib.loc = libname),"tcl", "tkScatterplotV3.tcl")
	tcl("source", tclfile)
	
}

.onAttach <- function(libname, pkgname) {
	packageStartupMessage("\nRnavGraph Version ",
			utils::packageDescription(pkg = pkgname, lib.loc = libname, field="Version"),
			'\nPlease read the package vignette. Use vignette("RnavGraph").\n\n')

        ## load Img tk extension if available
	sysname <- Sys.info()[1]
	didLoad <- TRUE
	if (sysname == "Darwin") {
            addTclPath("/System/Library/Tcl")
            didLoad <- tclRequire('Img')
	} else {
            didLoad <- tclRequire('Img')
	}
	
	if(identical(didLoad,FALSE)) {
            packageStartupMessage("Can not load the tk Img extension. Hence you can not use the 'ng_image_files' R function. Read the package vignette on how to set up tcl/tk.")	
	}
}
