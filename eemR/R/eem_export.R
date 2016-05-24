#' Export EEMs to Matlab
#'
#' @param file The .mat file name where to export the structure.
#' @param ... One or more object of class \code{eemlist}.
#'
#' @details The function exports EEMs into PARAFAC-ready Matlab \code{.mat} file
#'   usable by the \href{www.models.life.ku.dk/drEEM}{drEEM} toolbox.
#'
#' @return A structure named \code{OriginalData} is created and contains:
#'
#'   \describe{ \item{nSample}{The number of eems.} \item{nEx}{The number of
#'   excitation wavelengths.} \item{nEm}{The number of emission wavelengths.}
#'   \item{Ex}{A vector containing excitation wavelengths.} \item{Em}{A vector
#'   containing emission wavelengths.} \item{X}{A 3D matrix (nSample X nEx X
#'   nEm) containing EEMs.} }
#'
#'   \code{sample_name} The list of sample names (i.e. file names) of the
#'   imported EEMs.
#'
#' @export
#' @examples
#' file <- system.file("extdata/cary/", package = "eemR")
#' eem <- eem_read(file, recursive = TRUE)
#'
#' export_to <- paste(tempfile(), ".mat", sep = "")
#' eem_export_matlab(export_to, eem)

eem_export_matlab <- function(file, ...){

  eem <- list(...)

  list_classes <- unlist(lapply(eem, function(x) {class(x)}))

  stopifnot(all(list_classes %in% c("eem", "eemlist")),
            file.info(dirname(file))$isdir,
            grepl(".mat", basename(file)))

  eem <- eem_bind(...)

  ## Number of eem
  nSample <- length(eem)

  #---------------------------------------------------------------------
  # Check emission wavelengths
  #---------------------------------------------------------------------
  nEm <- unique(unlist(lapply(eem, function(x) length(x$em))))

  if(length(nEm) != 1){
    stop("Length of emission vectors are not the same across all eem.",
         call. = FALSE)
  }

  Em <- mapply(function(x) x$em, eem)

  if(ncol(unique(Em, MARGIN = 2)) != 1){
    stop("Emission vectors are not the same across all eem.",
         call. = FALSE)
  }

  Em <- Em[, 1] ## Just get the first column

  #---------------------------------------------------------------------
  # Check excitation wavelengths
  #---------------------------------------------------------------------
  nEx <- unique(unlist(lapply(eem, function(x) length(x$ex))))

  if(length(nEx) != 1){
    stop("Length of excitation vectors are not the same across all eem.",
         call. = FALSE)
  }

  Ex <- mapply(function(x) x$ex, eem)

  if(ncol(unique(Ex, MARGIN = 2)) != 1){
    stop("Exctiation vectors are not the same across all eem.",
         call. = FALSE)
  }

  Ex <- Ex[, 1] ## Just get the first column

  #---------------------------------------------------------------------
  # Prepare the 3D X matrix contianing eem sample nSample x nEm x nEx
  #---------------------------------------------------------------------

  ncol = unique(unlist(lapply(eem, function(x) ncol(x$x))))

  if(length(ncol) != 1){
    stop("EEMs do not have all the same number of columns across the dataset.",
         call. = FALSE)
  }

  nrow = unique(unlist(lapply(eem, function(x) nrow(x$x))))

  if(length(nrow) != 1){
    stop("EEMs do not have all the same number of rows across the dataset.",
         call. = FALSE)
  }

  X <- simplify2array(lapply(eem, function(x)x$x))

  X <- array(aperm(X, c(3, 1, 2)), dim = c(nSample, nEm, nEx))

  ## Use PARAFAC "naming" convention
  OriginalData <- list(X = X,
                       nEm = nEm,
                       nEx = nEx,
                       nSample = nSample,
                       Ex = Ex,
                       Em = Em)

  R.matlab::writeMat(file, OriginalData = OriginalData,
                     sample_names = eem_names(eem))

  message("Successfully exported ", nSample, " EEMs to ", file, ".\n")

}

