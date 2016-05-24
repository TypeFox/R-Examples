.get_co <- function(con) {
  seek(con, where = 0, origin = "start", rw = "r")
  readBin(con, what = character(), size = 3, endian = "little")
}

.get_comments <- function(con) {
  seek(con, where = 3, origin = "start", rw = "r")
  readBin(con, what = character(), size = 157, endian = "little")
}

.get_when <- function(con) {
  seek(con, 160)
  tm <- readBin(con, "integer", size = 2, n = 6, endian = "little")
  ISOdatetime(tm[6] + 1900, tm[5] + 1, tm[4], tm[3], tm[2], tm[1])
}

.get_program_version <- function(con) {
  seek(con, 178)
  res <- readBin(con, what = integer(), size = 1)
  major <- bitwShiftR(res, 4)
  minor <- bitwAnd(res, 7)
  paste0(major, '.', minor)
}

.get_file_version <- function(con) {
  seek(con, 179)
  res <- readBin(con, raw(), size = 1)
  res <- as.integer(res)
  major <- bitwShiftR(res, 4)
  minor <- bitwAnd(res, 7)
  paste0(major, '.', minor)
}

.get_dc_corr <- function(con) {
  seek(con, 181)
  readBin(con, logical(), size = 1)
}

.get_dc_time <- function(con) {
  seek(con, 182)
  res <- readBin(con, integer(), size = 4, endian = "little")
  as.POSIXct(res, origin = "1970-01-01")
}

.get_data_type <- function(con) {
  # from ASD doc
  types <- c('raw', 'reflectance', 'radiance', 'no_units', 'irradiance', 'qi', 'transmittance', 'unknown', 'absorbance')
  seek(con, 186)
  res <- readBin(con, integer(), size = 1)
  types[res + 1]
}

.get_ref_time <- function(con) {
  seek(con, 187)
  res <- readBin(con, integer(), size = 4, endian = "little")
  as.POSIXct(res, origin = "1970-01-01")
}

.get_ch1_wavel <- function(con) {
  seek(con, 191)
  readBin(con, numeric(), size = 4, endian = "little")
}

.get_wavel_step <- function(con) {
  seek(con, 195)
  readBin(con, numeric(), size = 4, endian = "little")
}

.get_data_format <- function(con) {
  # from ASD doc
  formats <- c('numeric', 'integer', 'double', 'unknown')
  seek(con, 199)
  res <- readBin(con, integer(), size = 1)
  formats[res - 1]
}

.get_channels <- function(con) {
  seek(con, 204)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_it <- function(con) {
  seek(con, 390)
  readBin(con, integer(), size = 4, endian = "little")
}

.get_fo <- function(con) {
  seek(con, 394)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_dcc <- function(con) {
  seek(con, 396)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_calibration <- function(con) {
  seek(con, 398)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_instrument_num <- function(con) {
  seek(con, 400)
  res <- readBin(con, integer(), size = 2, endian = "little")
  as.character(res)
}

.get_ip_numbits <- function(con) {
  seek(con, 418)
  readBin(con, integer(), size = 2, endian = "little")
}

# Flag codes:
# flags(0) AVGFIX'ed
# flags(1):
#   vnir saturation =1
#   swir1 satruation = 2
#   swir2 saturation = 3
#   Tec1 alarm= 8
#   Tec2 alarm = 16
#
.get_flags <- function(con) {
  seek(con, 421)
  readBin(con, integer(), size = 4, endian = "little")
}

.get_dc_count <- function(con) {
  seek(con, 425)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_ref_count <- function(con) {
  seek(con, 427)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_sample_count <- function(con) {
  seek(con, 429)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_instrument <- function(con) {
  instruments <- c('unknown', 'PSII', 'LSVNIR', 'FieldSpec VNIR', 'FieldSpec FR', 'FieldSpec NIR', 'CHEM', 'FieldSpec FullRange Unattended')
  seek(con, 431)
  res <- readBin(con, integer(), size = 1)
  instruments[res + 1]
}

.get_bulb <- function(con) {
  seek(con, 432)
  readBin(con, integer(), size = 4, endian = "little")
}

.get_swir1_gain <- function(con) {
  seek(con, 436)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_swir2_gain <- function(con) {
  seek(con, 438)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_swir1_offset <- function(con) {
  seek(con, 440)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_swir2_offset <- function(con) {
  seek(con, 442)
  readBin(con, integer(), size = 2, endian = "little")
}

.get_splice1_wavelength  <- function(con) {
  seek(con, 444)
  readBin(con, numeric(), size = 4, endian = "little")
}

.get_splice2_wavelength  <- function(con) {
  seek(con, 448)
  readBin(con, numeric(), size = 4, endian = "little")
}

.get_metadata <- function(con) {

  md <- list()

  # get file version
  md$co <- .get_co(con)
  # get comments
  md$comments <- .get_comments(con)
  # get time of acquisition
  md$when <- .get_when(con)
  # get version of program that created file
  md$program_version <- .get_program_version(con)
  # get spectrum file format version
  md$file_version <- .get_file_version(con)
  # flag dark correction
  md$dc_corr <- .get_dc_corr(con)
  # time of last dark correction
  md$dc_time <- .get_dc_time(con)
  # type of data stored
  md$data_type <- .get_data_type(con)
  # time of last white reference
  md$ref_time <- .get_ref_time(con)
  # starting wavelength in nm
  md$ch1_wavel <- .get_ch1_wavel(con)
  # wavelength step in nm
  md$wavel_step <- .get_wavel_step(con)
  # data format
  md$data_format <- .get_data_format(con)

  # ignoring the old_dc_count,old_ref_count,
  # old_sample count and application for now

  # number of channels in detector
  md$channels <- .get_channels(con)

  # ignoring the APP_DATA and GPS_DATA
  # for now

  # integration time, in ms
  md$it <- .get_it(con)
  # fore optic information
  md$fo <- .get_fo(con)
  # dark current correction value
  md$dcc <- .get_dcc(con)
  # calibration series
  md$calibration <- .get_calibration(con)
  # instrument number
  md$instrument_num <- .get_instrument_num(con)

  # ignoring the ymin, ymax, xmin, xmax for now

  # instrument dynamic range
  md$ip_numbits <- .get_ip_numbits(con)

  # ignoring xmode for now

  # warning flag
  md$flags <- .get_flags(con)
  # number of spectra averaged in the dark correction
  md$dc_count <- .get_dc_count(con)
  # number of spectra averaged in the white reference
  md$ref_count <- .get_ref_count(con)
  # number of spectra averaged in the created spectrum
  md$sample_count <- .get_sample_count(con)
  # instrument type that created the spectrum
  md$instrument <- .get_instrument(con)
  # id number of the bulb
  md$bulb <- .get_bulb(con)
  # gain for SWIR1
  md$swir1_gain <- .get_swir1_gain(con)
  # gain for SWIR2
  md$swir2_gain <- .get_swir2_gain(con)
  # offset for SWIR1
  md$swir1_offset <- .get_swir1_offset(con)
  # offset for SWIR2
  md$swir2_offset <- .get_swir2_offset(con)
  # wavelength of VNIR and SWIR1 splice
  md$splice1_wavelength <- .get_splice1_wavelength(con)
  # wavelength of SWIR1 and SWIR2 splice
  md$splice2_wavelength <- .get_splice2_wavelength(con)

  # ignoring SmartDetector data for now

  md
}

# constructs vector of wavelengths from metadata
.get_wavelengths <- function(md) {
  # compute the last wavelength
  start <- md$ch1_wavel
  step <- md$wavel_step
  n <- md$channels
  end <- start + n * step - 1
  # construct vector of wl
  seq(start, end, by = step)
}
