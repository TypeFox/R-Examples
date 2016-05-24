.onLoad <- function(libname, pkgname){
  #Set the options
  .setOptionsCurrent()
  .setOptionsDepreciated()
  
  #Set the theme and the last coordinates.
  theme_set(theme_gray())
}

#------------------------------------------------------------------------------
#CURRENT OPTIONS
#------------------------------------------------------------------------------
.setOptionsCurrent <- function(){
  options("tern.expand"                = 0.2)
  options('tern.margin'                = unit(0,'pt'))
  options("tern.default.T"             = "y")
  options("tern.default.L"             = "x")
  options("tern.default.R"             = "z")
  options("tern.clockwise"             = FALSE)
  options("tern.showtitles"            = TRUE)
  options("tern.showlabels"            = TRUE)
  options("tern.arrow.start"           = 0.3)
  options("tern.arrow.finish"          = 0.7)
  options('tern.arrowsep'              = 0.1)
  options('tern.vshift'                = 0.0)
  options('tern.hshift'                = 0.0)
  options("tern.showarrows"            = TRUE)
  options("tern.showgrid.major"        = TRUE)
  options("tern.showgrid.minor"        = TRUE)
  options("tern.ticks.outside"         = TRUE)
  options("tern.ticks.showprimary"     = TRUE)
  options("tern.ticks.showsecondary"   = FALSE)
  options("tern.breaks.default"        = seq(0.2, 1.0,by=0.2))
  options("tern.breaks.default.minor"  = seq(0.1, 0.9,by=0.2))
}

#------------------------------------------------------------------------------
#DEPRECIATED OPTIONS -- ie either not used anymore or in depreciated functions.
#------------------------------------------------------------------------------
.setOptionsDepreciated <- function(){
  options("tern.discard.external"      = TRUE)
  options("tern.expand.contour.inner"  =-0.0005)
  options("tern.dont_transform"        = FALSE)
  options("tern.mesh.buffer"           = 1.50)
  options("tern.mesh.size"             = 200)
}
