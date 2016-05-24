.onAttach <- function(...){
  ## figure out this year automatically 
  this.year <- substr(as.character(Sys.Date( )), 1, 4)

  packageStartupMessage("##\n## RxCEcolInf\n")
  packageStartupMessage("## Copyright (C) 2008-", this.year,
      " RxCEcolInf\n\n", sep="")

}

.onUnload <- function(libpath) {
    library.dynam.unload("RxCEcolInf", libpath)
}


#.First.lib <-
#    function (libname, pkgname)
#{
#    ## figure out this year automatically 
#    this.year <- substr(as.character(Sys.Date( )), 1, 4)
#    
#    ## echo output to screen
#    cat("##\n## RxCEcolInf\n")
#    cat("## Copyright (C) 2008-", this.year,
#        " RxCEcolInf\n\n", sep="")
#    require(MASS, quietly=TRUE)
#    require(mvtnorm, quietly=TRUE)
#    require(MCMCpack, quietly=TRUE)
#
#    library.dynam(pkgname, pkgname, lib.loc=libname)
#}






