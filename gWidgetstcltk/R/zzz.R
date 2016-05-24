.onLoad <- function(libname,pkgname,...) {
  ## methods isn't loaded yet, so we try calling through ::
  oldClasses <- c("tkwin", "tclVar", "tclObj")
  methods::setClass("tcltkObject")
  lapply(oldClasses, function(i) {
    methods::setOldClass(i)
    methods::setIs(i,"tcltkObject")
  })



  
}
         

       
tcltkStockIcons <- TcltkStockIcons$new()


.onAttach <- function(...) {
  ## version check
  if(as.numeric(.Tcl("info tclversion")) < 8.5) {
    packageStartupMessage("\n\n *** gWidgetstcltk needs tcl/tk version 8.5 or newer ***\n\n")
  }
  
  ## some configuration
  .Tcl("option add *tearOff 0")         # disable tearoff menus

  
  ## read in tklibs (from tcltk2 pacakge)
  f <- system.file("tklibs", "tablelist5.6", package="gWidgetstcltk")
  if(file.exists(f)) addTclPath(f)
  tclRequire("tablelist")
  sapply(c("tablelistConfig.tcl", "tablelistBind.tcl", "tablelistBind.tcl",
           "tablelistUtil.tcl", "tablelistEdit.tcl"), function(i) {
             f <-  system.file("tklibs", "tablelist5.6", "scripts", i, package="gWidgets2tcltk")
             if(file.exists(f))
               tcl("source", f)
          })

  f <- system.file("tklibs", "tooltip1.4", package="gWidgetstcltk")
  if(file.exists(f))
    addTclPath(f)
  try(tclRequire("tooltip"), silent=TRUE)
  f <- system.file("tklibs", "autoscroll.tcl", package="gWidgetstcltk")
  if(file.exists(f))
    tcl("source", f)
  


  ## ## read in tklibs (from tcltk2 pacakge)
  ## addTclPath(system.file("tklibs", package="gWidgetstcltk"))
  ## tclRequire("tooltip")
  ## tclRequire("autoscroll")


  ## Icons
  tcltkStockIcons$load_gWidgets_icons()
  ## use.table options
  ## images from http://ryanfait.com/resources/custom-checkboxes-and-radio-buttons/. Thanks
  f <- system.file("images", "checkbutton-off.gif", package="gWidgetstcltk")
  if(file.exists(f))
    tkimage.create("photo", "::image::off", file=f)
  f <- system.file("images", "checkbutton-on.gif",  package="gWidgetstcltk")
  if(file.exists(f))
    tkimage.create("photo", "::image::on",  file=f)
}
