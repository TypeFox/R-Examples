#' @import foreign Cairo RGtk2 VIM survey tkrplot gWidgetsRGtk2
#' @importFrom Hmisc spss.get
#' @importFrom Hmisc sasxport.get
.onAttach <- function(libname, pkgname) {
#    library.dynam("VIM", pkgname, libname)
    # load data
    #data(chorizonDL, package=pkgname, lib.loc=libname)
    #data(tao, package=pkgname, lib.loc=libname)
    #data(sleep, package=pkgname, lib.loc=libname)
    #data(kola.background, package=pkgname, lib.loc=libname)
    # prepare GUI
    tklibs <- file.path(find.package(pkgname, lib.loc=libname)[1], "tklibs")
    addTclPath(tklibs)
    tclRequire("BWidget")
    #vmGUImenu()
    packageStartupMessage("VIM is ready to use. \n (Enter vmGUImenu() to start the OLD graphical user interface.)\n
    Enter VIMGUI() to start the NEW graphical user interface.")
}

