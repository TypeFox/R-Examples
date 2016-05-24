# nmfgpu4R - R binding for the nmfgpu library
# 
# Copyright (C) 2015-2016  Sven Koitka (svenkoitka@fh-dortmund.de)
# 
# This file is part of nmfgpu4R.
# 
# nmfgpu4R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# nmfgpu4R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with nmfgpu4R.  If not, see <http://www.gnu.org/licenses/>.

#' R binding for computing non-negative matrix factorizations using CUDA
#' 
#' R binding for the libary \emph{nmfgpu} which can be used to compute Non-negative Matrix Factorizations (NMF) using CUDA hardware
#' acceleration. 
#' 
#' The main function to use is \code{\link{nmf}} which can be configured using various arguments.
#' In addition to it a few helper functions are provided, but they aren't necessary for using \code{\link{nmf}}.
#' 
#' @docType package
#' @name nmfgpu4R
#' @useDynLib nmfgpu4R
#' @import Rcpp
#' @import Matrix
#' @importFrom stats runif terms
#' @importFrom utils txtProgressBar unzip
NULL

nmfgpu4R.env <- new.env()
nmfgpu4R.env$initialized <- F

#' Initializes the C++ library nmfgpu, which provides the core functionality of this package.
#' 
#' @details As this package depends on a C++ library, there are some restrictions in terms of usage. First
#' this package is only compatible with x64 environments, because some CUDA libraries are not available for
#' x86 environments. Second you need a CUDA capable device starting with Kepler architecture, CUDA
#' device drivers and the CUDA toolkit with version 7.0 or higher. Lastly you need the nmfgpu library itself.
#' This package provides a basic service of downloading precompiled versions from github, if it is available for
#' your operating system and CUDA toolkit version. Otherwise you need to compile and install the library on
#' your own by following the instructions on \url{https://github.com/razorx89/nmfgpu}.
#' 
#' Even if the package downloads a precompiled version, it must not neccessarily be compatible with your system.
#' For example on Windows platforms you must have installed the "Microsoft Visual C++ Redistributable Packages for 
#' Visual Studio 2013", which can be found at \url{https://www.microsoft.com/en-us/download/details.aspx?id=40784}.
#' Furthermore on unix systems there could be a version mismatch with the \code{libstdc++.so} library, because the
#' installed compiler and the compiler which was used to build the binary could be different.
#' 
#' If you encounter any problems loading the nmfgpu library, then try to compile it by yourself.
#' 
#' @param quiet If true then informative messages about the found CUDA version and nmfgpu location will be suppressed.
#' 
#' @export
nmfgpu4R.init <- function(quiet=F) {
  # Check if it is already initialized
  if(nmfgpu4R.env$initialized) {
    warning("nmfgpu4R is already initialized, library will be reloaded")
    shutdownAdapters()
    nmfgpu4R.env$initialized <- F
  }
  
  # Check if it is a 64-bit R environment
  if(.Machine$sizeof.pointer != 8) {
    stop("nmfgpu4R is only compatible with 64-bit versions of R")
  }
  
  # Check if CUDA is installed
  cudaValid <- F
  cudaPath <- Sys.getenv("CUDA_PATH")
  if(cudaPath == "") {
    cudaPath <- Sys.getenv("CUDA_ROOT")
  }
  
  if(cudaPath != "") {
    # Try to find CUDA version
    output <- system2("nvcc", "--version", stdout=T)
    if(length(output) >= 4) {
      matches <- stringr::str_match(output[4], "^Cuda compilation tools, release ([0-9]+).([0-9]+)")
      if(all(!is.na(matches))) {
        cudaMajor <- as.numeric(matches[1,2])
        cudaMinor <- as.numeric(matches[1,3])
        cudaValid <- T
      }
    }
  }
  
  if(cudaValid) {
    if(!quiet) {
      message("CUDA toolkit:")
      message("  - Path: ", cudaPath)
      message("  - Version: ", cudaMajor, ".", cudaMinor)
    }
  } else {
    stop("No CUDA toolkit detected on this system")
  }
  
  # Check if CUDA is at least 7.0
  if(cudaMajor < 7) {
    stop("Library nmfgpu is only compatible with CUDA toolkit version 7.0 or higher")
  }
  
  # Check if nmfgpu has been installed by the user
  customBuild <- T
  nmfgpuRoot <- Sys.getenv("NMFGPU_ROOT")
  
  if(nmfgpuRoot == "") {
    nmfgpuRoot <- .downloadLibrary(cudaMajor, cudaMinor, quiet)
    if(is.null(nmfgpuRoot)) {
      return
    }
    customBuild <- F
  }
  
  # Output information about chosen location
  if(!quiet) {
    message("Loading nmfgpu from: ", nmfgpuRoot)
  }
  
  # 
  errmsg <- NULL
  if(!dir.exists(nmfgpuRoot)) {
    errmsg <- paste0("'", nmfgpuRoot, "' is not a valid and/or existing directory!")
  } else {
    if(!initializeAdapters(nmfgpuRoot)) {
      errmsg <- "Initialization of nmfgpu failed!"
      
      if(customBuild) {
        errmsg <- paste(errmsg, "Please recheck the installation instructions at https://github.com/razorx89/nmfgpu")
      } else {
        errmsg <- paste(errmsg, "As you are using a precompiled version of nmfgpu, please consider compiling it yourself. Maybe a dependency of the build toolchain is missing on your system.")
      }
    }
  }
  
  # Output error message if any is available
  if(!is.null(errmsg)) {
    stop(errmsg)
  }
  
  nmfgpu4R.env$initialized <- T
}

.ensureInitialized <- function() {
  if(!nmfgpu4R.env$initialized) {
    stop("Package is not initialized, please call nmfgpu4R.init() first!")
  }
}

# Unloads the C++ library when the package is unloaded
.onUnload <- function(lib, pkg) {
  if(nmfgpu4R.env$initialized) {
    shutdownAdapters()
  }
}

.getOperatingSystemIdentifier <- function() {
  if(.Platform$OS.type == "windows") { 
    "win"
  } else if(Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if(.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}

.downloadLibrary <- function(cudaMajor, cudaMinor, quiet) {
  # Generate URLs
  os <- .getOperatingSystemIdentifier()
  ver <- nmfgpuVersionString()
  libURL <- paste0("https://github.com/razorx89/nmfgpu/releases/download/", ver, "/nmfgpu_", os, "_cuda-", cudaMajor, "-", cudaMinor, ".zip")
  md5URL <- paste0(libURL, ".md5")
  
  # Generate filesystem paths
  pkgPath <- dirname(system.file(".", package = "nmfgpu4R"))
  installDir <- file.path(pkgPath, "inst", paste0("nmfgpu-", ver, "-cuda-v", cudaMajor, ".", cudaMinor))
  tmpArchive <- tempfile(fileext = ".zip")
  tmpMD5 <- tempfile(fileext = ".md5")
  
  # Check if directory exists
  if(dir.exists(installDir) && length(list.files(installDir, all.files=T, include.dirs=T, no..=T)) > 0) {
    return(installDir)
  } else {
    dir.create(path = installDir, showWarnings = F)
  }
  
  # Download from github
  if(!quiet) {
    message("Performing one-time download of library nmfgpu:")
    message("   URL: ", libURL)
  }
  tryCatch({
    utils::download.file(url = libURL, destfile = tmpArchive, mode = "wb", cacheOK = F, quiet = T)
    utils::download.file(url = md5URL, destfile = tmpMD5, mode = "w", cacheOK = F, quiet = T)
  }, error=function(cond) {
    warning("Precompiled version of nmfgpu ", ver, " for CUDA toolkit v", cudaMajor, ".", cudaMinor, " is not available. Please visit https://github.com/razorx89/nmfgpu and follow the installation instructions!")
    return(NULL)
  })
  
  if(!quiet) {
    message("Verifying integrity of archive file...")
  }
  md5Check <- tolower(readLines(tmpMD5, n=1, warn=F))
  md5Tmp <- tolower(as.character(tools::md5sum(tmpArchive)))
  
  if(md5Check != md5Tmp) {
    stop("Could not verify archive download, MD5 hashes do not match!")
  }
  
  # Extract archive
  unzip(tmpArchive, exdir=installDir)
  
  return(installDir)
}

# Sets a callback which gets called when a new frobenius norm has been calculated. Information like iteration number, frobenius norm, 
# RMSD and the delta frobenius norm are provided.
# 
# @param callbackFunction Either a valid R function to process data during algorithm execution or the NULL object to unset any callback.
#
# @note The callback function need to accept the following parameters: iteration, frobenius, deltaFrobenius, rmsd.
# 
# @export
#nmfgpu.setCallback <- function(callbackFunction) {
#  return(adapterSetCallback(callbackFunction))
#}


#' Requests the currently available and total amount of device memory.
#' 
#' @param deviceIndex If specified the memory info retrieval is restricted to the passed device indices. By default no restriction is active and
#' therefore memory information about all available CUDA devices are retrieved.
#' 
#' @return On success a list of lists will be returned, containing the following informations:
#' \tabular{ll}{
#'  \code{index} \tab Index of the CUDA device\cr
#'  \code{free.bytes} \tab Amount of free memory in bytes. \cr
#'  \code{total.bytes} \tab Total amount of memory in bytes.
#' }
#' 
#' @export
deviceMemoryInfo <- function(deviceIndex=NA) {
  .ensureInitialized()
  
  if(is.na(deviceIndex)) {
    deviceIndex = 0:(deviceCount()-1)
  }
  
  result <- list()
  for(i in deviceIndex) {
    tmp <- cppInfoForGpuIndex(i)
    tmp$index = i
    result <- c(result, list(tmp))
  }
  
  class(result) <- "DeviceMemoryInfo"
  
  return(result)
}

library(utils)

#' Prints the information of a 'DeviceMemoryInfo' object.
#' @param x Object of class 'DeviceMemoryInfo'
#' @param ... Other arguments
#' 
#' @export
#' @method print DeviceMemoryInfo
print.DeviceMemoryInfo <- function(x, ...) {
  for(i in 1:length(x)) {
    deviceInfo <- x[[i]]
    
    # Try to format bytes
    if(requireNamespace("gdata", quietly=T)) {
      used <- gdata::humanReadable(deviceInfo$total.bytes - deviceInfo$free.bytes, digits=2)
      total <- gdata::humanReadable(deviceInfo$total.bytes, digits=2)
    } else {
      used <- paste(deviceInfo$total.bytes - deviceInfo$free.bytes, "B")
      total <- paste(deviceInfo$total.bytes, "B")
    }
    
    cat("#", deviceInfo$index, ": ", deviceInfo$name, " [ Allocated: ", used, " / ", total, " ]\n", sep="")
    pb <- txtProgressBar(initial=(deviceInfo$total.bytes - deviceInfo$free.bytes) / deviceInfo$total.bytes, style=3)
  }
}

#' Retrieves the total number of installed CUDA devices. 
#' @export
deviceCount <- function() {
  .ensureInitialized()
  return(cppNumberOfGpu())
}

#' Selects the specified device as primary computation device. All further invocations to nmfgpu will use the specified
#' CUDA device. 
#' 
#' @param deviceIndex Index of the CUDA device, which should be used for computation.
#' 
#' @note CUDA enumerates devices starting with 0 for the first device. 
#'  
#' @export
chooseDevice <- function(deviceIndex) {
  .ensureInitialized()
  
  if(!is.numeric(deviceIndex) || deviceIndex %% 1 != 0 || deviceIndex < 0) {
    stop("device.index must be a non-negative integer number")
  }
  
  if(!cppChooseGpu(deviceIndex)) {
    stop("Failed to choose specified device!")
  }
}