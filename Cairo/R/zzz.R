.onLoad <- function(libname, pkgname) {
    ## add our libs to the PATH
    if (.Platform$OS.type=="windows") {
        lp<-gsub("/","\\\\",paste(libname,pkgname,"libs",sep="/"))
        cp<-strsplit(Sys.getenv("PATH"),";")
        if (! lp %in% cp) Sys.setenv(PATH=paste(lp,Sys.getenv("PATH"),sep=";"))
    }
    library.dynam("Cairo", pkgname, libname)
    .Call("Rcairo_initialize", PACKAGE="Cairo")
}
