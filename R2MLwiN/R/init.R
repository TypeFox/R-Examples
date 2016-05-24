.onAttach <- function(...) {
  if (is.null(getOption("MLwiN_path"))) {
    if (Sys.info()["sysname"] == "Linux") {
      options(MLwiN_path = "/usr/bin/mlnscript")
    }
    if (Sys.info()["sysname"] == "Darwin") {
      options(MLwiN_path = "/opt/mln/mlnscript")
    }
    if (Sys.info()["sysname"] == "FreeBSD") {
      options(MLwiN_path = "/usr/local/bin/mlnscript")
    }
    if (Sys.info()["sysname"] == "Windows") {
      options(MLwiN_path = "C:/Program Files (x86)/MLwiN v2.36/")
    }
  }
  packageStartupMessage(paste0("The MLwiN_path option is currently set to ", getOption("MLwiN_path"), "\n", "To change this use: options(MLwiN_path=\"<path to MLwiN>\")\n"))
} 
