getSysinfo=function()
##title<< Compile system information
##description<< getSysinfo is a convenience function to compile some information about the
## computing system and environment used.
##details<< The function is mainly used to save the system environment information
##          in ncdf files containing the results of some calculations. 
{
    package.names=sapply(sessionInfo()[['otherPkgs']],'[[','Package')
    package.versions=sapply(sessionInfo()[['otherPkgs']],'[[','Version')
    packages.all=paste(package.names,package.versions,collapse=', ',sep='')
    pars.sys=c('user','nodename','sysname','release')
    R.system=paste(sessionInfo()[[1]]$version.string)
    sys.info = paste(pars.sys,Sys.info()[pars.sys],collapse=', ',sep=':')
    all.info=paste(c(sys.info,', ',R.system,', installed Packages: ',packages.all),sep='',collapse='')
    ##value<< character string with all version and system information of the current R system
    return(all.info)
}
