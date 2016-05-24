#' @include get-spectrum.r

#' @title Reads ASD Binary Files in R.
#'
#' @description \code{asdreader} implements a simple reader to read spectroscopy data collected using ASD (now PAN Analytics, Inc.) visible near-infrared spectrometers, and stored using the ASD format (which is documented here: \url{http://support.asdi.com/Document/Viewer.aspx?id=95}.
#'
#' The spectra can be extracted from the ASD file as raw (DN), white reference, radiance, or reflectance. Additionally, the metadata information contained in the ASD file header can also be accessed.
#'
#' @docType package
#' @author Pierre Roudier
#' @name asdreader
NULL

#' @name get_metadata
#' @export
#' @title Reads metadata header from ASD file
#' @author Pierre Roudier
#' @param f character, path to ASD file
#' @return a list storing the metadata information in the ASD header,
#' as documented here: \url{http://support.asdi.com/Document/Viewer.aspx?id=95}
#' @examples
#' asd_fn <- asd_file()
#' md <- get_metadata(asd_fn)
#' names(md)
#'
get_metadata <- function(f) {
  # Open connection to file
  con <- file(f, "rb")

  # Get metadata from file
  md <- .get_metadata(con)

  # Close connection
  close(con)

  # Return list
  md
}

#' @name get_spectra
#' @title Reads reflectance from ASD file
#' @export
#' @author Pierre Roudier
#' @param f a vector of paths to ASD file(s)
#' @param type a character vector, which type of spectra to return. \code{"reflectance"}, \code{"raw"}, \code{"white_reference"} are currently supported
#' @return a matrix of the spectrum contained in the ASD file(s)
#' @examples
#'
#' # Get the path to the demo file
#'
#' asd_fn <- asd_file()
#' print(asd_fn)
#'
#' # Example with one file name
#'
#' m1 <- get_spectra(asd_fn)
#' matplot(t(m1), type = 'l')
#'
#' # Example with a vector of file names
#'
#' asd_fns <- rep(asd_fn, times = 4)
#' print(asd_fns) # (in this case, 4 times the same file)
#'
#' m2 <- get_spectra(asd_fns)
#' matplot(t(m2), type = 'l')
#'
get_spectra <- function(f, type = "reflectance") {

  if (length(type) > 1) {
    stop('Please select only one type.')
  }

  res <- lapply(f, function(fn) {

    # Open connection to file
    con <- file(fn, "rb")

    # Get metadata from file
    md <- .get_metadata(con)

    # Read spectrum
    spec <- .get_spec(con, md)

    # Close connection
    close(con)

    # Process spectra depending on data type required
    res <- .process_spectra(spec, md, type)

    # Return spectrum with proper column names
    res <- matrix(res, ncol = length(res))
    colnames(res) <- .get_wavelengths(md)

    res
  })

  do.call(rbind, res)
}
