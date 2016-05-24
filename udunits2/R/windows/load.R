.onLoad<-function(libname,pkgname){
    Sys.setenv(UDUNITS2_XML_PATH = file.path(libname, pkgname, "share/udunits/udunits2.xml"))
}
