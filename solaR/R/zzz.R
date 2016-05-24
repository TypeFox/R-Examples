.solEnv <- new.env()
assign('oldTZ', Sys.getenv('TZ'), envir = .solEnv)

.onLoad <- function(libpath, pkgname)
{
    Sys.setenv(TZ='UTC')
}

.onAttach <- function(libpath, pkgname){
    packageStartupMessage('Time Zone set to UTC.\n')
}
    
.onUnload <- function(libpath)
{
    oldTZ <- get('oldTZ', envir = .solEnv)
    Sys.setenv(TZ = oldTZ)
}
