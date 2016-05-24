## Adapted from Rcpp

## Make sure system.file returns an absolute path
Rclusterpp.system.file <- function(...){
    tools::file_path_as_absolute( base::system.file( ..., package = "Rclusterpp" ) )
}

## Identifies if the default linking on the platform should be static
## or dynamic. Currently only linux uses dynamic linking by default
## although it works fine on mac osx as well
staticLinking <- function() {
    ! grepl( "^linux", R.version$os )
}

## Use R's internal knowledge of path settings to find the lib/ directory
## plus optinally an arch-specific directory on system building multi-arch
RclusterppLdPath <- function() {
    if (nzchar(.Platform$r_arch)) {	## eg amd64, ia64, mips
        path <- Rclusterpp.system.file("lib",.Platform$r_arch)
    } else {
        path <- Rclusterpp.system.file("lib")
    }
    path
}

## Provide linker flags -- i.e. -L/path/to/libRclusterpp -- as well as an
## optional rpath call needed to tell the Linux dynamic linker about the
## location.  This is not needed on OS X where we encode this as library
## built time (see src/Makevars) or Windows where we use a static library
## Updated Jan 2010:  We now default to static linking but allow the use
##                    of rpath on Linux if static==FALSE has been chosen
##                    Note that this is probably being called from LdFlags()
RclusterppLdFlags <- function(static=staticLinking()) {
    rclusterppdir <- RclusterppLdPath()
    if (static) {                               # static is default on Windows and OS X
        flags <- paste(rclusterppdir, "/libRclusterpp.a", sep="")
    } else {					# else for dynamic linking
        flags <- paste("-L", rclusterppdir, " -lRclusterpp", sep="") # baseline setting
        if ((.Platform$OS.type == "unix") &&    # on Linux, we can use rpath to encode path
            (length(grep("^linux",R.version$os)))) {
            flags <- paste(flags, " -Wl,-rpath,", rclusterppdir, sep="")
        }
    }
    invisible(flags)
    #invisible(paste(Rcpp::RcppLdFlags(static), flags, sep=" "))
}

## LdFlags defaults to static linking on the non-Linux platforms Windows and OS X
LdFlags <- function(static=staticLinking()) {
    cat(RclusterppLdFlags(static=static))
}
