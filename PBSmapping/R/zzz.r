# Taking cue from Roger Bivand's maptools:
.PBSmapEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onLoad <- function(lib, pkg)
{
	library.dynam("PBSmapping", pkg, lib);
}

.onAttach <- function(lib, pkg)
{
	assign("PBSprint",FALSE,envir=.PBSmapEnv)

	# obtain values necessary for the start-up message
	pkg_info <- utils::sessionInfo( package="PBSmapping" )$otherPkgs$PBSmapping
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	userguide_path <- system.file( "doc/PBSmapping-UG.pdf", package = "PBSmapping")
	year <- substring(date(),nchar(date())-3,nchar(date()))
	
	packageStartupMessage("
-----------------------------------------------------------
PBS Mapping ", pkg_info$Version, " -- Copyright (C) 2003-",year," Fisheries and Oceans Canada

PBS Mapping comes with ABSOLUTELY NO WARRANTY;
for details see the file COPYING.
This is free software, and you are welcome to redistribute
it under certain conditions, as outlined in the above file.

A complete user guide 'PBSmapping-UG.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
http://code.google.com/p/pbs-software/

To see demos, type '.PBSfigs()'.
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSmapEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	"bcBathymetry",
	"deldir",
	"nepacLL","nepacLLhigh",
	"PBSval","pythagoras",
	"read.dbf",
	"surveyData",
	"towData","towTracks",
	"worldLL"),
	package="PBSmapping")

