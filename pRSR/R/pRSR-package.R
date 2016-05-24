".onLoad" <-
function(lib, pkg){ library.dynam("pRSR", pkg, lib)
data("TPout", package=pkg, envir=parent.env(environment())) 
}
