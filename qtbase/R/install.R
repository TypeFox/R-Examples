.windows_qt_path <- function()
  file.path(system.file(package = "qtbase"), "qt", .Platform$r_arch)

.install_system_dependencies <- function()
{
  windows32_config <-
    list(
         source = FALSE,
         url =
         "http://ggobi-documentation.googlecode.com/files/qt-4.8.0-win32.zip",
         installer = function(path) {
           dest_path <- .windows_qt_path()
           ## unzip does this, but we want to see any warnings
           dir.create(dest_path, recursive = TRUE) 
           unzip(path, exdir = dest_path)
         }
         )

  windows64_config <- NULL
  
  darwin_config <- NULL
  
  unix_config <- NULL
  
  web <- "http://qt.nokia.com"
  
  install_system_dep <- function(dep_name, dep_url, dep_web, installer)
    {
      if (!interactive()) {
        message("Please install ", dep_name, " from ", dep_url)
        return()
      }
      choice <- menu(paste(c("Install", "Do not install"), dep_name), TRUE, 
                     paste("The qtbase package requires", dep_name))
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
           "Please install Qt manually, if necessary. See: ", web)
    
    install_system_dep("Qt", config$url, web, config$installer)
  }
  
  install_all()
}
