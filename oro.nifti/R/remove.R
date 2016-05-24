#' @title Remove File Extensions Around the NIfTI/ANALYZE Formats
#' 
#' @description Simple function(s) that remove file extensions commonly found when using
#' NIfTI-1 or ANALYZE format files.
#' 
#' 
#' @aliases rmniigz rmnii rmgz rmhdrgz rmhdr rmimggz rmimg
#' @param x is the file name.
#' @return The file name without offending suffix.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @keywords Misc
#' @name rmniigz
#' @rdname remove
#' @export
rmniigz <- function(x) {
  sub(".nii.gz", "", x)
}
#' @rdname remove
#' @export
rmnii <- function(x) {
  sub(".nii", "", x)
}
#' @rdname remove
#' @export
rmgz <- function(x) {
  sub(".gz", "", x)
}
#' @rdname remove
#' @export
rmhdrgz <- function(x) {
  sub(".hdr.gz", "", x)
}
#' @rdname remove
#' @export
rmhdr <- function(x) {
  sub(".hdr", "", x)
}
#' @rdname remove
#' @export
rmimggz <- function(x) {
  sub(".img.gz", "", x)
}
#' @rdname remove
#' @export
rmimg <- function(x) {
  sub(".img", "", x)
}

