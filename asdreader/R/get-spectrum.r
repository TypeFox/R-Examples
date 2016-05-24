#' @include get-metadata.r

.normalise_spectrum <- function(spec, md) {

  wl <- .get_wavelengths(md)

  idx1 <- which(wl <= md$splice1_wavelength)
  idx2 <- which(wl > md$splice1_wavelength & wl <= md$splice2_wavelength)
  idx3 <- which(wl > md$splice2_wavelength)

  spec[idx1] <- spec[idx1] / md$it
  spec[idx2] <- spec[idx2] * md$swir1_gain / 2048
  spec[idx3] <- spec[idx3] * md$swir2_gain / 2048

  spec
}

# get spectra
.get_spec <- function(con, md) {
  seek(con, 484)
  # Raw spectrum
  spec <- readBin(con, what = md$data_format, n = md$channels, endian = "little")

  # White reference flag
  # seek(con, 17692) = 484 + 8 * md$channels
  wr_flag <- readBin(con, what = logical(), size = 1)

  # White reference time
  # seek(con, 17693) = 484 + 8 * md$channels + 1
  wr_time <- readBin(con, integer(), size = 8, endian = "little")

  # Spectrum time
  # seek(con, 17701) = 484 + 8 * md$channels + 9
  spec_time <- readBin(con, integer(), size = 8, endian = "little")

  # Spectrum description length
  # seek(con, 17709) = 484 + 8 * md$channels + 17
  spec_description_length <- readBin(con, integer(), size = 2, endian = "little")
  # Spectrum description
  # seek(con, 17710) # = 484 + 8 * md$channels + 19
  spec_description <- readBin(con, character(), size = spec_description_length, endian = "little")

  # White reference
  # seek(con, 17712) = 484 + 8 * md$channels + 20
  wr <- readBin(con, what = md$data_format, n = md$channels, endian = "little")

  res <- list(spectrum = spec, wr = wr)
}

.process_spectra <- function(spec, md, type) {
  if (type == 'reflectance') {
    if (md$data_type == 'reflectance') {
      res <- spec$spectrum
    } else if (md$data_type == 'raw') {
      res <- .normalise_spectrum(spec$spectrum, md) / .normalise_spectrum(spec$wr, md)
    } else {
      stop(paste0('File only contains data of type ', md$data_type, '.'))
    }
  } else if (type == 'radiance') {
    if (md$data_type == 'radiance') {
      res <- spec$spectrum
    } else if (md$data_type == 'raw') {
      res <- .normalise_spectrum(spec$spectrum, md)
    } else {
      stop(paste0('File only contains data of type ', md$data_type, '.'))
    }
  } else if (type == 'raw') {
    if (md$data_type == 'raw') {
      res <- spec$spectrum
    } else {
      stop(paste0('File only contains data of type ', md$data_type, '.'))
    }
  } else if (type == 'white_reference') {
    res <- .normalise_spectrum(spec$wr, md)
  } else {
    stop('Invalid type.')
  }
}

