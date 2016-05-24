#' Function to install the Blat software
#'
#' @param arch Architechture. Either 32 or 64.
#' @param force logical. Force installation or ask for confirmation first?
install_blat <- function(arch = 64, force = FALSE)
{  
  if (!arch %in% c(32, 64))
    stop("Use either 32 or 64 as value for arch argument.", call. = FALSE)
  
  install_msg <- sprintf(
    "This will install the %dbit version of Blat in %s.\nContinue? [Yes/Y/No/N]",
    arch, system.file(package = .packageName))
  
  if (!force) {
    message(install_msg)
    if (!tolower(readLines(n = 1L)) %in% c("y", "ye", "yes")) {
      message("Nothing is installed. Exiting.")
      return(invisible(NULL))
    }
  }
  
  url <- c(
    `32` = "https://github.com/smbache/blatr/raw/master/blat323_32.full.zip", 
    `64` = "https://github.com/smbache/blatr/raw/master/blat323_64.full.zip"
    )
  
  temp_dir   <- tempdir()
  temp_file  <- file.path(temp_dir, "blat.zip")
  dest_dir   <- system.file(package = .packageName)
  blat_files <- paste0("blat323/full/blat", c(".dll", ".exe", ".lib", "dll.h"))
  
  if (download.file(url[as.character(arch)], 
                    temp_file, quiet = TRUE) != 0)
    stop("Could not download the files. Try manual installation.", call. = FALSE)
  
  unzip(temp_file, blat_files, exdir = dest_dir, junkpaths = TRUE)   
    
  if (!check_install())
    stop("The installation failed.", call. = FALSE)
  
  message("The installation was successful.")
  
  invisible(NULL)
}

#' Get the needed Blat files.
#'
#' @return a character vector of needed Blat files.
blat_files <- function() 
{
  system.file(c("blat.exe", "blat.dll", "blat.lib", "blatdll.h"), 
              package = .packageName)
}

#' Function to check whether Blat is installed correctly.
#'
#' @return logical.
check_install <- function()
{
  !any(blat_files() == "")
}
