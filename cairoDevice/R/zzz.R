.onLoad <- function(libname, pkgname)
{
  if (.Platform$OS.type == "windows") {
    dllpath <- file.path(.windows_gtk_path(), "bin")
    dll <- try(library.dynam("cairoDevice", pkgname, libname,
                             DLLpath = dllpath),
               silent = getOption("verbose"))
  }
  else dll <- try(library.dynam("cairoDevice", pkgname, libname),
                  silent = getOption("verbose"))
  
  if (is.character(dll)) {
    warning("Failed to load cairoDevice, attempting to install it", dll)
    .install_system_dependencies()
    return()
  }
  
  #library.dynam("cairoDevice", pkgname, libname)
  .C("loadGTK", success = logical(1), PACKAGE="cairoDevice")

  # register device as being interactive
  deviceIsInteractive("Cairo")
}

.closeDevices <- function()
{
    devices <- dev.list()
    gtk.devices <- devices[names(devices)=="Cairo"]
    if(length(gtk.devices) > 0) {
        dev.off(gtk.devices)
    }
}

.onUnload <- function(libpath) {
  .closeDevices()
  .C("cleanupGTK", PACKAGE = "cairoDevice")
  library.dynam.unload("cairoDevice", libpath)
}

.windows_gtk_path <- function()
  file.path(system.file(package = "cairoDevice"), "gtk", .Platform$r_arch)

.install_system_dependencies <- function()
{
  windows32_config <-
    list(
         source = FALSE,
         gtk_url = "http://ftp.gnome.org/pub/gnome/binaries/win32/gtk+/2.22/gtk+-bundle_2.22.1-20101227_win32.zip",
         installer = function(path) {
           gtk_path <- .windows_gtk_path()
           ## unzip does this, but we want to see any warnings
           dir.create(gtk_path, recursive = TRUE) 
           utils::unzip(path, exdir = gtk_path)
         }
         )

  windows64_config <- windows32_config
  windows64_config$gtk_url <- "http://ftp.gnome.org/pub/gnome/binaries/win64/gtk+/2.22/gtk+-bundle_2.22.1-20101229_win64.zip"
  
  darwin_config <- list(
                        source = FALSE,
                        gtk_url = "http://r.research.att.com/libs/GTK_2.18.5-X11.pkg", 
                        installer = function(path) {
                          system(paste("open", path))
                        }
                        )
  
  unix_config <- NULL
  
  gtk_web <- "http://www.gtk.org"
  
  install_system_dep <- function(dep_name, dep_url, dep_web, installer)
    {
      if (!interactive()) {
        message("Please install ", dep_name, " from ", dep_url)
        return()
      }
      choice <- utils::menu(paste(c("Install", "Do not install"), dep_name), T, 
                            paste("Need", dep_name,
                                  "? (Restart R after installing)"))
      if (choice == 1) {
        path <- file.path(tempdir(), basename(sub("\\?.*", "", dep_url)))
        if (utils::download.file(dep_url, path, mode="wb") > 0)
          stop("Failed to download ", dep_name)
        installer(path)
      }
      message("Learn more about ", dep_name, " at ", dep_web)
    }
  
  install_all <- function() {
    if (.Platform$OS.type == "windows") {
      if (.Platform$r_arch == "i386")
        config <- windows32_config
      else config <- windows64_config
    } else if (length(grep("darwin", R.version$platform))) 
      config <- darwin_config
    else config <- unix_config
    
    if (is.null(config))
      stop("This platform is not yet supported by the automatic installer. ",
           "Please install GTK+ manually, if necessary. See: ", gtk_web)
    
    install_system_dep("GTK+", config$gtk_url, gtk_web, config$installer)
  }
  
  install_all()
  
  message("If the package still does not load, please ensure that GTK+ is",
          " installed and that it is on your PATH environment variable")
  message("IN ANY CASE, RESTART R BEFORE TRYING TO LOAD THE PACKAGE AGAIN")
}
