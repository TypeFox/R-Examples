if (is.R()){

    ".onLoad" <- function(lib, pkg){
        if (.Platform$r_arch == "i386") {
            .onLoad.i386(lib, pkg)
        }
        else if (.Platform$r_arch == "x64"){
            .onLoad.x64(lib, pkg)
        }
        else {
            stop("Unknown architecture ", .Platform$r_arch, " , should be i386 or x64")
        }
    }

    ".onUnload" <- function(libpath){
        if (.Platform$r_arch == "i386") {
            .onUnload.i386(libpath)
        }
        else if (.Platform$r_arch == "x64"){
            .onUnload.x64(libpath)
        }
        else {
            stop("Unknown architecture ", .Platform$r_arch, " , should be i386 or x64")
        }
    }    

}
