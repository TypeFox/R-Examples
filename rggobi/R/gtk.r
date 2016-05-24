# RGtk casting
# Cast an object to an RGtk object
# 
# @keyword internal 
as.RGtkObject <- function(x) {
  class(x) <- c(class(x), "RGtkObject")
  x
}

# Gtk main window
# Retrieve RGtk object for main window
#
# Useful for embedding in other applications or for listening
# to their signals via RGtk2.
# 
# @keyword internal 
ggobi_gtk_main_window <- function(.gobi = ggobi_get()) {
  .GGobiCall("getMainWindow", .gobi = .gobi)
}


# Gtk menu bar
# Retrieve RGtk object for menu bar
# 
# Useful for embedding in other applications or for listening
# to their signals via RGtk2.
# 
# @keyword internal 
ggobi_gtk_menu_bar <- function(.gobi = ggobi_get()) {
  .GGobiCall("getMenubar", .gobi = .gobi)
}

