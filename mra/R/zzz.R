.onAttach<-function(libname, pkgname){

	#library.dynam("mra", pkgname, libname)  # deprecated
    #v <- packageDescription("mra", fields="Version")  # this requires utils package, and I don't want to make mra dependent on utils, just for this.

    packageStartupMessage( "Mark-Recapture Analysis (vers 2.16.4)" )  # You have to change this every version
	#packageStartupMessage("\ncontact author: Trent McDonald (tmcdonald@west-inc.com)\ncontributions: Eric Regehr, Jeff Bromaghin") 
}


