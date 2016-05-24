#' Get a polgon giving area of coverage of a given WRS-1 or WRS-2 path and row
#'
#' @export
#' @import wrspathrowData
#' @param wrs_path WRS-1 or WRS-2 path as an integer
#' @param wrs_row WRS-1 or WRS-2 row as an integer
#' @param wrs_type 1 (for WRS-1) or 2 (for WRS-2)
#' @param wrs_mode either 'D' for descending (daytime) or 'A' for ascending 
#' (nighttime)
#' @return list with path and row as integers
#' @examples
#' library(sp)
#'
#' x <- pathrow_poly(225, 61)
#' plot(x)
pathrow_poly <- function(wrs_path, wrs_row, wrs_type='2', wrs_mode='D') {
    if (wrs_type == 2) {
        if (wrs_path < 1 || wrs_path > 233) {
            stop('WRS-2 paths range from 1 to 233')
        }
        if (wrs_row < 1 || wrs_row > 248) {
            stop('WRS-2 rows range from 1 to 248')
        }
    } else if (wrs_type == 1) {
        if (wrs_path < 1 || wrs_path > 251) {
            stop('WRS-1 paths range from 1 to 251')
        }
        if (wrs_row < 1 || wrs_row > 248) {
            stop('WRS-1 rows range from 1 to 248')
        }
    }
    wrs_polys <- load_wrs_data(wrs_type, wrs_mode)
    return(wrs_polys[wrs_polys@data$PATH == wrs_path &
                     wrs_polys@data$ROW == wrs_row, ])
}

