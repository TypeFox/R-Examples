.gtkArgs <-
function()
{
 c("R")
}

.onLoad <-
function(libname, pkgname)
{
 options(depwarn = TRUE, gdkFlush = TRUE)

 if (.Platform$OS.type == "windows") {
   dllpath <- Sys.getenv("RGTK2_GTK2_PATH")
   if (!nzchar(dllpath))
     dllpath <- file.path(.windows_gtk_path(), "bin")
   dll <- try(library.dynam("RGtk2", pkgname, libname, DLLpath = dllpath),
              silent = getOption("verbose"))
 }
 else dll <- try(library.dynam("RGtk2", pkgname, libname),
                 silent = getOption("verbose"))
 if (is.character(dll)) {
   warning("Failed to load RGtk2 dynamic library, attempting to install it.",
           call. = FALSE)
   .install_system_dependencies()
   if (.Platform$OS.type == "windows") # just try to load the package again
     .onLoad(libname, pkgname)
   return()
 }
   

 if(is.function(.gtkArgs))
  args <- .gtkArgs()
 else
  args <- as.character(.gtkArgs)

 if(!(gtkInit(args))) {
   .init_failed()
 }

 .initClasses()
}

.onUnload <- function(libpath) {
  .gtkCleanup()
  library.dynam.unload("RGtk2", libpath)
}

.windows_gtk_path <- function()
  file.path(system.file(package = "RGtk2"), "gtk", .Platform$r_arch)

.configure_gtk_theme <- function(theme) {
  ## Only applies to Windows so far
  config_path <- file.path(system.file(package = "RGtk2"), "gtk",
                           .Platform$r_arch, "etc", "gtk-2.0")
  dir.create(config_path, recursive = TRUE)
  writeLines(sprintf("gtk-theme-name = \"%s\"", theme),
             file.path(config_path, "gtkrc"))
}

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
           unzip(path, exdir = gtk_path)
           .configure_gtk_theme("MS-Windows")
         }
         )

  windows64_config <- windows32_config
  windows64_config$gtk_url <- "http://ftp.gnome.org/pub/gnome/binaries/win64/gtk+/2.22/gtk+-bundle_2.22.1-20101229_win64.zip"
  
  darwin_config <- list(
    source = FALSE,
    gtk_url = "http://r.research.att.com/libs/GTK_2.24.17-X11.pkg", 
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
    choice <- menu(paste(c("Install", "Do not install"), dep_name), T, 
      paste("Need", dep_name, "? (Restart R after installing)"))
    if (choice == 1) {
      path <- file.path(tempdir(), basename(sub("\\?.*", "", dep_url)))
      if (download.file(dep_url, path, mode="wb") > 0)
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

.init_failed <- function() {
  message("R session is headless; GTK+ not initialized.")
}
