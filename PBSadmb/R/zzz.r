# Taking cue from Roger Bivand's maptools:
.PBSadmbEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

#.onLoad <- function(libname, pkgname)
.onAttach <- function(libname, pkgname)
{
	.initOptions()
	pkg_info <- utils::sessionInfo( package="PBSadmb" )$otherPkgs$PBSadmb
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	
	userguide_path <- system.file( "doc/PBSadmb-UG.pdf", package = "PBSadmb" )
	year <- substring(date(),nchar(date())-3,nchar(date()))
	
	packageStartupMessage("
-----------------------------------------------------------
PBS ADMB ", pkg_info$Version, " -- Copyright (C) 2008-",year," Fisheries and Oceans Canada

A complete user guide 'PBSadmb-UG.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
http://code.google.com/p/pbs-software/

Type admb() to start a GUI for operating ADMB.
-----------------------------------------------------------

")
}

.onUnload <- function(libpath) {
	rm(.PBSadmbEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	".PBSadmb",".PBSadmb.pkgOptions",".PBSmod",
	"add","admbpath","admbver","argvec",
	"chkadmb","chkgcc",
	"debugsymbols","digest","dll",
	"editor",
	"files",
	"gccpath","gccver",
	"isdir",
	"logfile",
	"nsims","nthin",
	"optfile",
	"PBSadmb","pltView","prefix","pthin",
	"raneff","runType",
	"safe",
	"toView",
	"verbose"
	), package="PBSadmb")

