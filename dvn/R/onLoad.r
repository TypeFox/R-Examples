.onLoad <- function(libname, pkgname){
	options(dvn = 'https://thedata.harvard.edu/dvn/') # set default dataverse
    options(dvn.user = '') # initialize dataverse username option
    options(dvn.pwd = '') # initialize dataverse password option
}