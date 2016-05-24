.PBSmodEnv <- new.env(FALSE, parent=globalenv())  # Taking cue from Roger Bivand's maptools

.onLoad <- function(libname, pkgname)
{
	library.dynam("PBSmodelling", pkgname, libname)
}

.onAttach <- function(libname, pkgname){
        # obtain values necessary for the start-up message
	pkg_info <- utils::sessionInfo( package="PBSmodelling" )$otherPkgs$PBSmodelling
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()
	userguide_path <- system.file( "doc/PBSmodelling-UG.pdf", package = "PBSmodelling" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Modelling ", pkg_info$Version, " -- Copyright (C) 2005-",year," Fisheries and Oceans Canada

A complete user guide 'PBSmodelling-UG.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
http://code.google.com/p/pbs-software/
-----------------------------------------------------------

")
	#Load custom PBSmodelling tcl scripts
	tcl("lappend", "auto_path", system.file( "tcl_scripts", package = "PBSmodelling" ) )
	tclRequire( "PBSmodelling" )

        #Ensure that Tk is available, first
        tk <- tclRequire("Tk", warn = FALSE)
        if ( is.logical (tk) ) {
		packageStartupMessage("
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERROR: Your system failed to load Tk, the widget toolkit required for
PBSmodelling.  See the R FAQ for your operating system at
        http://cran.r-project.org/faqs.html
for suggestions on resolving issues with Tcl/Tk.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        }
        
        #Try to load the required packages; set warn = FALSE so that
        #the user won't see warnings during the install; warnings
        #during the install don't include the the details (below), and
        #thus, are very confusing
	bwidget <- tclRequire("BWidget", warn = FALSE)
        if ( is.logical (bwidget) ) {
                # try the included distribution
		tcl("lappend", "auto_path", system.file( "thirdparty/BWidget-1.9.0/", package = "PBSmodelling" ) )
		bwidget <- tclRequire("BWidget", warn = FALSE)
        }
	tktable <- tclRequire("Tktable", warn = FALSE)

        #If Tk succeeded and either package failed to load, display an appropriate error message
        if( !is.logical (tk) && ( is.logical( bwidget ) || is.logical( tktable ) ) ) {
		packageStartupMessage("
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERROR: Your system is missing at least one Tcl/Tk package that
PBSmodelling requires to function correctly (details below).  Missing
packages must be installed outside of R.")
                if( is.logical( bwidget ) ) {
                        packageStartupMessage("
MISSING PACKAGE: \"BWidget\" (Tcl/Tk package)

Ubuntu/Debian users may try installing BWidget with the command
        sudo apt-get install bwidget

If this command doesn't work, the source code for the BWidget
toolkit can be downloaded from their SourceForge site
        http://sourceforge.net/projects/tcllib/files/BWidget/")
                }
                if( is.logical( tktable ) ) {
                        packageStartupMessage("
MISSING PACKAGE: \"Tktable\" (Tcl/Tk package)

Ubuntu/Debian users may try installing Tktable with the command
        sudo apt-get install libtktable2.10
or
        sudo apt-get install tk-table
If these commands fail, it may be possible to obtain a version of the
package \"tk-table\" from
        http://packages.ubuntu.com
and then install it with a command like
        sudo dpkg -i tk-table_2.10-1_i386.deb
If the above (dpkg) command fails because the package \"tk\" cannot be
found, first try installing the package \"tk\" with the command
        sudo apt-get install tk

Mac users may try installing the Tcl/Tk package available from
        http://cran.r-project.org/bin/macosx/tools
This package includes Tktable 2.9.  Those Mac users who are using
MacPorts (http://www.macports.org/) may try installing Tktable with
the command
        sudo port install tktable

If these commands don't work, the source code for tktable can be
downloaded from their SourceForge site
        http://sourceforge.net/projects/tktable/files/tktable/")
                }

		packageStartupMessage("
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	}
        
	.initPBSoptions()
}
.onUnload <- function(libpath) {
	rm(.PBSmodEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	".cfl",".closeChoice",".cls",".cwd",".dls",".dwd",".makeChoice",".PBSmod",".runExHelperQuit",
	"blist",
	"chFile","chTest","closeALL","closeSDE","command","curVal","cwd",
	"for",
	"global",
	"httpdPort",
	"maxVal","minVal",
	"OK",
	"PBS.history","PBSmin","prefix",
	"remote","runs",
	"tgot","tmp.before",
	"v.tab","variable",
	"wN"
	), package="PBSmodelling" )

