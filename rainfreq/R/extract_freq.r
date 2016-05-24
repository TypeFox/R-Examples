#' @title Extract rainfall frequency estimates for desired region of the USA
#'
#' @param region_name region name. choose one of se(southeast), sw(southwest), 
#' mw(midwest), orb(ohio river basin and surrounding areas), hi(hawaii) 
#' and ak(alaska). default is se(southeast). \pkg{rainfreq} regional selection 
#' criterion is currently limited to the 50 states (plus DC). 
#' @param storm_RP storm return period in years. choose one of 1, 2, 5, 10, 25,
#' 50, 100, 200, 500, 1000. default is 100 years.
#' @param storm_duration storm duration in minutes, hours, or days. choose one
#' of 5m, 10m, 15m, 30m, 60m (in minutes); 2h, 3h, 6h, 12h, 24h, 48h 
#' (in hours); 3d, 4d, 7d, 10d, 30d, 45d, 60d (in days). default is 24h.
#' @param flag_down_read flag to indicate both downloading and reading of the 
#' data is desired. default is TRUE.
#' @param flag_down_only flag to indicate only downloading of the data is 
#' desired. default is FALSE.
#' @param flag_read_only flag to indicate only reading of the data is 
#' desired. default is FALSE. if set to TRUE a valid file path is required.
#' @param nws_data_path location of downloaded nws zip files. when 
#' flag_read_only is TRUE it defaults to the working directory
#' 
#' @return RasterLayer, if flag_down_only is set to FALSE, NULL otherwise; if 
#' the NWS website is not working a value of 10 is returned
#' 
#' @author Gopi Goteti
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # southeast, 100yr-24hour storm
#' x_se <- extract_freq()
#' class(x_se)
#' print(x_se)
#' # midwest, 1000yr-48hour storm
#' x_mw <- extract_freq(region_name = "mw", storm_RP = 1000, storm_duration = "48h")
#' print(y_se)
#' # download only, southeast, 100yr-24hour storm
#' extract_freq(flag_down_read = FALSE, flag_down_only = TRUE)
#' # read after download, southeast, 100yr-24hour storm
#' x_se <- extract_freq(flag_down_read = FALSE, flag_read_only = TRUE)
#' print(x_se)
#' }
#' # record rainfall for the usa
#' data(rain_max_usa)
#' head(rain_max_usa)
#' # record rainfall for the world
#' data(rain_max_world)
#' head(rain_max_world)
extract_freq <- function(region_name = "se", 
                         storm_RP = 100, 
                         storm_duration = "24h", 
                         flag_down_read = TRUE,
                         flag_down_only = FALSE,
                         flag_read_only = FALSE,
                         nws_data_path = "") {
  
  # check for region names
  name_vec <- c("sw", "orb", "mw", "se", "ak", "hi")
  if (!(region_name %in% name_vec)) {
    stop(paste("region_name: region name must be one of - ", paste(name_vec, collapse = ",")))
  }
  
  # check return periods
  rp_vec <- c(1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
  if (!(storm_RP %in% rp_vec)) {
    stop(paste("storm_RP: return period must be one of - ", paste(rp_vec, collapse = ",")))
  }
  
  # check duration
  duration_vec <- c("5m", "10m", "15m", "30m", "60m", "2h", "3h", "6h", "12h", "24h", "48h", "3d", "4d", "7d", "10d", "20d", "30d", "45d", "60d")
  if (!(storm_duration %in% duration_vec)) {
    stop(paste("storm_duration: duration must be one of - ", paste(duration_vec, collapse = ",")))
  }
  
  # check data access flags
  if (!is.logical(flag_down_read)) {
    stop("flag_down_read should be either TRUE or FALSE!")
  }  
  if (!is.logical(flag_down_only)) {
    stop("flag_down_only should be either TRUE or FALSE!")
  }
  if (!is.logical(flag_read_only)) {
    stop("flag_read_only should be either TRUE or FALSE!")
  }
  # only one flag could be true
  if (sum(flag_down_read, flag_down_only, flag_read_only) != 1) {
    stop("one and only one of flag_down_read, flag_down_only, flag_read_only could be TRUE!")
  }

  # construct url and file name
  duration1 <- substr(storm_duration, 1, nchar(storm_duration) - 1)
  duration1 <- sprintf("%.2d", as.integer(duration1))
  duration2 <- substr(storm_duration, nchar(storm_duration), nchar(storm_duration))
  new_duration <- paste0(duration1, duration2)  
  
  file_dump <- paste0(region_name, storm_RP, "yr", new_duration, "a.zip")
  
  if (flag_down_read | flag_down_only) {
    
    # NWS ftp site
    data_source <- "ftp://hdsc.nws.noaa.gov/pub/hdsc/data/"
    
    # check whether the site exists/responds
    if (!url.exists(data_source)) {
      message(paste("NWS website", data_source, "is not responding, please try later!"))
      return (10)
    }
    url_name <- paste0(data_source, region_name, "/", file_dump)
    if (!url.exists(url_name)) {
      message(paste("NWS website", url_name, "is not responding or does not exist, 
                 please check and/or try later!"))
      return (10)
    }
    
    # download using RCurl
    message("downloading data. might take a few moments ...")
    file_dump_con <- CFILE(file_dump, "wb")
    curlPerform(url = url_name, writedata = file_dump_con@ref, verbose = FALSE)
    close(file_dump_con)
  }  
  
  if (flag_down_read | flag_read_only) {
    
    # ascii file within the zipped file has a different suffix
    if (region_name %in% c("sw", "se", "ak", "mw")) {
      file_asc <- gsub(".zip", ".asc", file_dump)  
    } else {
      file_asc <- gsub(".zip", ".ver3", file_dump)  
    }
    
    if (flag_read_only) {
      # append file location path
      if (nws_data_path == "") {
        nws_data_path <- getwd()
      }
      file_dump <- file.path(nws_data_path, file_dump)
      if (!file.exists(file_dump)) {
        stop(paste("file", file_dump, "does not exist in", nws_data_path))
      }
    }
    
    # open connection to the ascii file inside the zipped file    
    message("reading contents of the zipped file ...")
    file_asc_con <- unz(description = file_dump, filename = file_asc)
    
    # extract header info, below code based on SDMTools read.asc function
    hdr_list <- scan(file_asc_con, what = list('', ''), nlines = 6, quiet = TRUE)
    hdr_info <- hdr_list[[2]]
    
    # number of columns, check for decimals
    if(grepl("\\.", hdr_info[1])) {
      stop("num cols must be an integer! cannot be a real number!")
    } else {
      nc <- as.integer(hdr_info[1])
    }
    
    # number of rows, check for decimals
    if(grepl("\\.", hdr_info[2])) {
      stop("num rows must be an integer! cannot be a real number!")
    } else {
      nl <- as.integer(hdr_info[2])
    }
    
    if(any(c(nc, nl) <= 0)) {
      stop("number of columns and rows must both be positive integers!")
    }
    
    # lower left corner, x
    xll <- as.double(hdr_info[3])
    # lower left corner, y
    yll <- as.double(hdr_info[4])
    
    # cell size
    cs <- as.double(hdr_info[5])
    if (cs <= 0) {
      stop("cell size/resolution should be a positive number!")
    }
    
    # nodata value
    if(grepl("\\.", hdr_info[6])) {
      nodata <- as.double(hdr_info[6])
    } else {
      nodata <- as.integer(hdr_info[6])
    }
    
    # read remainder of the file    
    output <- scan(file_asc_con, nmax = nl * nc, skip = 6, quiet = TRUE)
    close(file_asc_con)
    
    #convert no data to NA
    output[output == nodata] <- NA
    
    # convert data to matrix, and re-order or flip as needed
    output <- matrix(c(as.matrix(output)), ncol = nl, nrow = nc)
    output <- output[, ncol(output):1]
    
    # define output's class attribs, consistent with SDMTools and adehabitat
    attr(output, "xll") <- xll
    attr(output, "yll") <- yll
    attr(output, "cellsize") <- cs
    attr(output, "type") <- 'numeric'
    class(output) <- "asc"
    
    # convert to raster
    output <- raster.from.asc(output)
    
    output
  }
}
