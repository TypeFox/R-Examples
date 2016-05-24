#' @importFrom RCurl getURL
#' @importFrom stringr str_extract
verify_download <- function(tile_url, local_path) {
    header <- getURL(tile_url, nobody=1L, header=1L)
    header <- strsplit(header, "\r\n")[[1]]
    content_length <- header[grepl('Content-Length: ', header)]
    remote_size <- as.numeric(str_extract(content_length, '[0-9]+'))
    local_size <- file.info(local_path)$size
    if (remote_size != local_size) {
        return(3)
    } else {
        return(0)
    }
}

#' @importFrom utils download.file
download_tile <- function(tile_url, local_path) {
    ret_code <- download.file(tile_url, local_path, mode="wb")
    if (ret_code != 0) {
        message(paste('Warning: problem downloading', basename(local_path)))
        return(1)
    } else if (verify_download(tile_url, local_path) != 0) {
        message(paste("Warning: verification failed on", basename(local_path)))
        return(2)
    } else {
        return(0)
    }
}

#' Download a set of GFC tiles
#'
#' This function first checks whether each tile in a set GFC product tiles is 
#' present locally, and that local file sizes match the file sizes of the files 
#' available on the Google server hosting the GFC product. Next, the function 
#' downloads all tiles that either are not present locally, or that are present 
#' but have file sizes differing from the file on the Google server.
#'
#' @seealso \code{\link{extract_gfc}}
#'
#' @export
#' @importFrom sp bbox
#' @importFrom stringr str_extract
#' @importFrom utils file_test
#' @param tiles \code{SpatialPolygonsDataFrame} with GFC 
#' product tiles to download, as calculated by the \code{calc_gfc_tiles} 
#' function.
#' @param output_folder the folder to save output data in
#' @param images which images to download. Can be any of 'treecover2000', 
#' 'loss', 'gain', 'lossyear', 'datamask', 'first', and 'last'.
#' @param data_year which version of the Hansen data to use
#' @examples
#' \dontrun{
#' output_folder <- 'H:/Data/TEAM/GFC_Product'
#' tiles <- calc_gfc_tiles(test_poly)
#' download_tiles(tiles, output_folder)
#' }
download_tiles <- function(tiles, output_folder,
                           images=c('treecover2000', 'loss', 'gain', 
                                    'lossyear', 'datamask'),
                           data_year=2015) {
    stopifnot(all(images %in% c('treecover2000', 'loss', 'gain', 'lossyear', 
                                'datamask', 'first', 'last')))
    if (!file_test('-d', output_folder)) {
        stop('output_folder does not exist')
    }
    message(paste(length(tiles), 'tiles to download/check.'))
    successes <- 0
    failures <- 0
    skips <- 0

    for (n in 1:length(tiles)) {
        gfc_tile <- tiles[n,]
        min_x <- bbox(gfc_tile)[1, 1]
        max_y <- bbox(gfc_tile)[2, 2]
        if (min_x < 0) {
            min_x <- paste0(sprintf('%03i', abs(min_x)), 'W')
        } else {
            min_x <- paste0(sprintf('%03i', min_x), 'E')
        }
        if (max_y < 0) {
            max_y <- paste0(sprintf('%02i', abs(max_y)), 'S')
        } else {
            max_y <- paste0(sprintf('%02i', max_y), 'N')
        }
        file_root <- paste0('Hansen_GFC', data_year, '_')
        file_suffix <- paste0('_', max_y, '_', min_x, '.tif')
        filenames <- paste0(file_root, images, file_suffix)

        tile_urls <- paste0(paste0('http://commondatastorage.googleapis.com/earthenginepartners-hansen/GFC', data_year, '/'), filenames)
        local_paths <- file.path(output_folder, filenames)

        for (i in 1:length(filenames)) {
            tile_url <- tile_urls[i]
            local_path <- local_paths[i]
            if (file.exists(local_path)) {
                if (verify_download(tile_url, local_path)) {
                    message(paste(basename(local_path), "exists but doesn't match remote - re-downloading file"))
                } else {
                    #message(paste(basename(local_path), 'exists and matches remote - skipping download'))
                    skips <- skips + 1
                    next
                }
            }
            if (download_tile(tile_url, local_path) == 0) {
                successes <- successes + 1
            } else {
                failures <- failures + 1
            }

        }
    }

    message(paste(successes, "file(s) succeeded,", skips, "file(s) skipped,", 
                  failures, "file(s) failed."))
}
