#.First.lib <- function(libname, pkgname) {
# library.dynam("TUWmodel", pkgname, libname)
#}
.onLoad <- function(libname, pkgname){
 # do whatever needs to be done when the package is loaded
 # some people use it to bombard users with 
 # messages using 
 #packageStartupMessage( "my package is so cool" )
 #packageStartupMessage( "so I will print these lines each time you load it")
 library.dynam("TUWmodel", pkgname, libname)
}

