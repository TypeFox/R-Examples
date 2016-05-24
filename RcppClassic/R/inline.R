
inlineCxxPlugin <- Rcpp:::Rcpp.plugin.maker( 
    include.before = "#include <RcppClassic.h>", 
    package = "RcppClassic", 
    Makevars = NULL, 
    Makevars.win = NULL, 
    libs = RcppClassicLdFlags()
    )
    
RcppClassic.system.file <- function(...) {
    tools::file_path_as_absolute(base::system.file(..., package = "RcppClassic"))
}

RcppClassicLdPath <- function() {
    if (nzchar(.Platform$r_arch)) {	## eg amd64, ia64, mips
        path <- RcppClassic.system.file("lib",.Platform$r_arch)
    } else {
        path <- RcppClassic.system.file("lib")
    }
    path
}

RcppClassicLdFlags <- function(static=staticLinking()) {
    rcppclassicdir <- RcppClassicLdPath()
    if (static) {                               # static is default on Windows and OS X
        flags <- paste(rcppclassicdir, "/libRcppClassic.a", sep="")
    } else {					# else for dynamic linking
        flags <- paste("-L", rcppclassicdir, " -lRcppClassic", sep="") # baseline setting
        if ((.Platform$OS.type == "unix") &&    # on Linux, we can use rpath to encode path
            (length(grep("^linux",R.version$os)))) {
            flags <- paste(flags, " -Wl,-rpath,", rcppclassicdir, sep="")
        }
    }
    invisible(flags)
}

LdFlags <- function(static=staticLinking()) {
    cat(RcppClassicLdFlags(static=static))
}

staticLinking <- function() {
    ! grepl( "^linux", R.version$os )
}

