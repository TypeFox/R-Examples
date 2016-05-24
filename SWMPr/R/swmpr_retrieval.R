######
#' Import current station records from the CDMO
#' 
#' Import current station records from the CDMO starting with the most current date
#' 
#' @param station_code chr string of station, 7 or 8 characters
#' @param Max numeric value for number of records to obtain from the current date
#' @param param chr string for a single parameter to return, defaults to all parameters for a station type.
#' @param trace logical indicating if import progress is printed in console
#' 
#' @export
#' 
#' @import httr
#' 
#' @seealso \code{\link{all_params_dtrng}}, \code{\link{single_param}}
#' 
#' @concept retrieve
#' 
#' @return  Returns a swmpr object, all available parameters including QAQC columns
#' 
#' @details 
#' This function retrieves data from the CDMO through the web services URL.  The computer making the request must have a registered IP address.  Visit the CDMO web services page for more information: \url{http://cdmo.baruch.sc.edu/webservices.cfm}.  Function is the CDMO equivalent of \code{exportAllParamsXMLNew} but actually uses \code{\link{all_params_dtrng}}, which is a direct call to \code{exportAllParamsDateRangeXMLNew}.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ## all parameters for a station, most recent
#' all_params('hudscwq')
#' 
#' }
all_params <- function(station_code, Max = 100, param = NULL, trace = TRUE){

  # url
  serv <- "http://cdmo.baruch.sc.edu/webservices2/requests.cfc?wsdl"
  
  # get from most recent record
  dat <- try({
    httr::GET(serv, 
      query = list(
        method = 'exportAllParamsXMLNew',
        station_code = station_code, 
        recs = 1
      )
    )
  }, silent = TRUE)
  
  # stop if retrieval error
  if('try-error' %in% class(dat))
    stop('Error retrieving data, check metadata for station availability.')

  # parse reply from server 
  dat <- parser(dat)
  
  # starting date as character
  dtrng <- dat$datetimestamp
  dtrng <- strsplit(as.character(dtrng), ' ')[[length(dtrng)]][1]
  dtrng <- c('01/01/1970', dtrng)
    
  # pass to all_params_dtrng
  out <- all_params_dtrng(station_code, dtrng, param = param, trace = trace, Max = Max)

  return(out)
  
}

######
#' Get CDMO records within a date range
#' 
#' Get station records from the CDMO within a date range
#' 
#' @param station_code chr string of station, 7 or 8 characters
#' @param dtrng two element chr string, each of format MM/DD/YYYY
#' @param param chr string for a single parameter to return, defaults to all parameters for a station type.
#' @param trace logical indicating if import progress is printed in console
#' @param Max numeric indicating maximum number of records to return
#' 
#' @export
#' 
#' @concept retrieve
#' 
#' @return Returns a swmpr object, all parameters for a station type (nutrients, water quality, or meteorological) or a single parameter if specified.  QAQC columns are not provided for single parameters.
#' 
#' @details 
#' This function retrieves data from the CDMO through the web services URL.  The computer making the request must have a registered IP address.  Visit the CDMO web services page for more information: \url{http://cdmo.baruch.sc.edu/webservices.cfm}.  This function is the CDMO equivalent of \code{exportAllParamsDateRangeXMLNew}.
#' Download time may be excessive for large requests.
#' 
#' @seealso \code{\link{all_params}}, \code{\link{single_param}}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ## get all parameters within a date range
#' all_params_dtrng('apaebwq', c('01/01/2013', '02/01/2013'))
#' 
#' ## get single parameter within a date range
#' all_params_dtrng('apaebwq', c('01/01/2013', '02/01/2013'), 
#'    param = 'do_mgl')
#' 
#' }
all_params_dtrng <- function(station_code, dtrng, param = NULL, trace = TRUE, Max = NULL){
  
  ##
  # access CDMO web services

  # timer
  tictoc::tic()
  
  # url
  serv <- "http://cdmo.baruch.sc.edu/webservices2/requests.cfc?wsdl"

  # for initializing while loop
  tz <- time_vec(station_code = station_code, tz_only = TRUE)
  obsmin <- as.POSIXct(dtrng[1], format = '%m/%d/%Y', tz = tz)
  end_obs <- obsmin + 1
  
  # object to fill
  out_all <- NULL
  
  if(trace) cat('Importing...\n\n')
  
  # start loop
  while(end_obs > obsmin){
    
    # arguments to pass to function on server
    web_args = list(
        method = 'exportAllParamsDateRangeXMLNew',
        station_code = station_code,
        mindate = dtrng[1],
        maxdate = dtrng[2]
        )
    
    # add a parameter argument if provided
    if(!is.null(param)) web_args$param <- param
    
    # request data
    dat <- try({
      httr::GET(
        serv,
        query = web_args
      )
    }, silent = TRUE)
    
    # stop if retrieval error
    if('try-error' %in% class(dat))
      stop('Error retrieving data, check metadata for station availability.')
    
    # parse reply from server 
    out <- parser(dat)

    # sometimes data request is good, but empty data frame returned
    if(nrow(out) == 0)
      stop('Empty data frame, check metadata for station availability')
    
    # type of parameter requested - wq, nut, or met, NOT the param argument
    parm <- substring(station_code, 6)
    nms <- param_names(parm)[[parm]]
    
    # format datetimestamp, sort, get relevant columns as data frame
    out[, 'datetimestamp'] <- time_vec(out[, 'datetimestamp'], station_code)
    out <- out[order(out$datetimestamp), ]
    out <- data.frame(
      datetimestamp = out$datetimestamp,
      out[, tolower(names(out)) %in% nms, drop = FALSE],
      row.names = 1:nrow(out)
      )
    names(out) <- tolower(names(out))

    # get new loop 
    end_obs <- min(out$datetimestamp)
    max_obs <- max(out$datetimestamp)

    # progress
    if(trace) cat('\t', as.character(as.Date(max_obs) - 1), 'to', as.character(as.Date(end_obs)), '\n')
    
    # exit if no new data
    if(!is.null(out_all)){
      if(end_obs == min(out_all$datetimestamp)) break
    }

    # append to output
    out_all <- rbind(out_all, out)
    out_all <- unique(out_all[order(out_all$datetimestamp), ])

    if(!is.null(Max)){
      if(nrow(out_all) >= Max){  
        out_all <- out_all[(1 + nrow(out_all) - Max):nrow(out_all), ]
        break
      }
    }
    
    # new date ranges
    dtrng[2] <- as.character(as.Date(end_obs))
    dtrng[2] <- paste0(substr(dtrng[2], 6, nchar(dtrng[2])), '/', substr(dtrng[2], 1, 4))
    dtrng[2] <- gsub('-', '/', dtrng[2])
    
  }
  
  # sort by date, then remove duplicates (there will be overlaps)
  # data columns as numeric
  parms <- nms[!grepl('^f_', nms)]
  out <- out_all
  out[, names(out) %in% parms] <- apply(out[, names(out) %in% parms, drop = FALSE], 2, as.numeric)
  row.names(out) <- 1:nrow(out)
  
  # convert to swmpr class
  out <- swmpr(out, station_code)
  
  if(trace){
    cat('\n')
    cat(nrow(out), 'records, ')
    tictoc::toc()
    }
  
  # return output
  return(out)
  
}

######
#' Get CDMO records for a single parameter
#' 
#' Get stations records from the CDMO for a single parameter starting with the most current date
#' 
#' @param station_code chr string of station, 7 or 8 characters
#' @param Max numeric value for number of records to obtain from the current date
#' @param param chr string for a single parameter to return.
#' @param trace logical indicating if import progress is printed in console
#' 
#' @import XML
#' 
#' @concept retrieve
#' 
#' @export
#' 
#' @return Returns a swmpr object with one parameter.  QAQC columns are not provided.
#' 
#' @details 
#' This function retrieves data from the CDMO through the web services URL. The computer making the request must have a registered IP address.  Visit the CDMO web services page for more information: \url{http://cdmo.baruch.sc.edu/webservices.cfm}.  This function is the CDMO equivalent of \code{exportSingleParamXML}.
#' 
#' @seealso \code{\link{all_params}}, \code{\link{all_params_dtrng}}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' ## single parameter for a station, most recent
#' single_param('hudscwq', 'do_mgl')
#' 
#' }
single_param <- function(station_code, param, Max = 100, trace = TRUE){

  ##
  # access CDMO web services
  
  # url
  serv <- "http://cdmo.baruch.sc.edu/webservices2/requests.cfc?wsdl"
  
  # request data
  dat <- try({
    httr::GET(
      serv,
      query = list(
        method = 'exportSingleParamXMLNew',
        station_code = station_code,
        param = param,
        recs = 1
      )
    )
  }, silent = TRUE)
  
  # stop if retrieval error
  if('try-error' %in% class(dat))
    stop('Error retrieving data, check metadata for station availability.')
  
  # parse reply from server 
  dat <- parser(dat)
  
  # sometimes data request is good, but empty data frame returned
  if(nrow(dat) == 0)
    stop('Empty data frame, check metadata for station availability')
  
  # starting date as character
  dtrng <- dat$datetimestamp
  dtrng <- strsplit(as.character(dtrng), ' ')[[length(dtrng)]][1]
  dtrng <- c('01/01/1970', dtrng)
    
  # pass to all_params_dtrng
  out <- all_params_dtrng(station_code, dtrng, param = param, trace = trace, Max = Max)

  return(out)
  
  # convert to swmpr class
  out <- swmpr(out, station_code)
  
  # return output
  return(out)
  
}

######
#' Import local CDMO data
#' 
#' Import local data that were obtained from the CDMO through the zip downloads feature
#' 
#' @param  path chr string of full path to .csv files with raw data, can be a zipped or unzipped directory where the former must include the .zip extension
#' @param  station_code chr string of station to import, typically 7 or 8 characters including wq, nut, or met extensions, may include full name with year, excluding file extension
#' @param  trace logical indicating if progress is sent to console, default \code{FALSE}
#' 
#' @concept retrieve
#' 
#' @export
#' 
#' @importFrom utils read.csv unzip
#' 
#' @return Returns a swmpr object with all parameters and QAQC columns for the station.  The full date range in the raw data are also imported.
#' 
#' @details 
#' The function is designed to import local data that were downloaded from the CDMO outside of R. This approach works best for larger data requests, specifically those from the zip downloads feature in the advanced query section of the CDMO. The function may also work using data from the data export system, but this feature has not been extensively tested. The downloaded data will be in a compressed folder that includes multiple .csv files by year for a given data type (e.g., apacpwq2002.csv, apacpwq2003.csv, apacpnut2002.csv, etc.). The import_local function can be used to import files directly from the compressed folder or after the folder is decompressed.  In the former case, the requested files are extracted to a temporary directory and then deleted after they are loaded into the current session.  An example dataset is available online to illustrate the format of the data provided through the zip downloads feature.  See the link below to access these data.  All example datasets included with the package were derived from these raw data.
#' 
#' Occasionally, duplicate time stamps are present in the raw data.  The function handles duplicate entries differently depending on the data type (water quality,  weather, or nutrients).  For water quality and nutrient data, duplicate time stamps are simply removed.  Note that nutrient data often contain replicate samples with similar but not duplicated time stamps within a few minutes of each other.  Replicates with unique time stamps are not removed but can be further processed using \code{\link{rem_reps}}.  Weather data prior to 2007 may contain duplicate time stamps at frequencies for 60 (hourly) and 144 (daily) averages, in addition to 15 minute frequencies.  Duplicate values that correspond to the smallest value in the frequency column (15 minutes) are retained.  
#' 
#' Zip download request through CDMO: \url{http://cdmo.baruch.sc.edu/aqs/zips.cfm}
#' 
#' Example dataset: \url{https://s3.amazonaws.com/swmpexdata/zip_ex.zip}
#' 
#' @seealso \code{\link{all_params}}, \code{\link{all_params_dtrng}}, \code{\link{rem_reps}}, \code{\link{single_param}}
#' 
#' @examples
#' 
#' \dontrun{
#' ## this is the path for csv example files, decompressed
#' path <- 'C:/this/is/my/data/path'
#'
#' ## import, do not include file extension
#' import_local(path, 'apaebmet') 
#' 
#' ## this is the path for csv example files, zipped folder
#' path <- 'C:/this/is/my/data/path.zip'
#'
#' ## import, do not include file extension
#' import_local(path, 'apaebmet') 
#' }
import_local <- function(path, station_code, trace = FALSE){
  
  # add .zip if not present
  if(file.exists(paste0(path, '.zip'))){
      path <- paste0(path, '.zip')
    }
  
  # check if file exists 
  if(!file.exists(path)){
    stop('Path does not exist')
  }
  
  # check if path is zipped
  zips <- grepl('\\.zip$', path)
  
  # remove file extension if present, lower case
  station_code <- tolower(gsub('\\.csv$', '', station_code))
  
  ##
  # find station files in path
  
  # for zipped
  if(zips){
    
    # get the file names in the zipped folder
    # check if the requested files exist
    file_nms <- unzip(path, list = TRUE)$Name
    expr <- paste0(station_code, '.*', '\\.csv$')
    files_in <- grep(expr, file_nms, value = TRUE, ignore.case = TRUE)
    if(length(files_in) == 0) stop('File(s) not found.')
    
    # extract to temporary file
    tmp_fl <- tempfile()
    unzip(path, files = files_in, exdir = tmp_fl)
    files_in <- dir(tmp_fl, recursive = TRUE)
    
    # reassign path to temporary file
    path <- tmp_fl
    
  # for unzipped    
  } else {

    file_nms <- dir(path)
    expr <- paste0('^', station_code, '.*', '\\.csv$')
  
  }

  files_in <- grep(expr, file_nms, value = TRUE, ignore.case = TRUE)
  if(length(files_in) == 0) stop('File(s) not found.')

  station_code <- tolower(station_code)
  
  # import all data files for a station
  dat <- vector('list', length(files_in))
  names(dat) <- gsub('.csv', '', files_in)
  
  if(trace) cat('Loading files...\n\n')
  
  for(file_in in files_in){
    
    if(trace) cat(file_in, '\t')
    
    ##
    # import
 
    # import file, try using read.csv, else readlines
    tmp <- try({
      read.csv(file.path(path, file_in), stringsAsFactors = FALSE)
    }, silent = TRUE)
    
    if('try-error' %in% class(tmp)){
      raw <- readLines(file.path(path, file_in))
      keep_lines <- grep(paste0('^', station_code), raw)
      tmp <- raw[keep_lines]
      tmp <- strsplit(tmp, ',')
      tmp <- do.call('rbind', tmp)
      tmp <- data.frame(tmp, stringsAsFactors = FALSE)
      names(tmp)  <- strsplit(
        gsub('["\\"]', '', raw[keep_lines[1] - 1]),
        ',')[[1]] 
    }
    
    names(tmp) <- tolower(names(tmp))
    
    # remove stationcode, isswmp columns
    tmp <- tmp[, !names(tmp) %in% c('stationcode', 'isswmp')]
    
    # convert date time to posix
    names(tmp)[grep('datetimestamp', names(tmp), ignore.case = TRUE)] <- 'datetimestamp'
    tmp$datetimestamp <- time_vec(tmp$datetimestamp, station_code)
    
    # append to output list
    nm <-  gsub('.csv', '', file_in)
    dat[[nm]] <- tmp
    
    }
  
  # remove temporary files if zips
  if(zips) unlink(tmp_fl, recursive = TRUE)
  
  ##
  # column names for each parameter type, used to subset combined data
  # kept as upper case here because improted data will match, changed to lower below

  # names to use
  parm <- substring(station_code, 6)
  parm <- gsub('[0-9.*]', '', parm)
  nms <- param_names(parm)[[parm]]
  
  ##
  # deal with duplicate time stamps depending on data type
  
  out <- do.call('rbind', dat)
  
  # if duplicated timestamps and met, keep those with minimum value in frequency
  if('met' %in% parm & any(duplicated(out$datetimestamp)) & 'frequency' %in% names(out)){
    
    min_step <- as.character(min(as.numeric(unique(out$frequency))))
    out <- out[out$frequency %in% min_step, ]
    
    # sometimes duplicates still remain at same frequency
    out <- out[!duplicated(out$datetimestamp),]  
    
  }
  
  # remove duplicate time stamps from wq and nut data
  if(any(c('nut', 'wq') %in% parm) & any(duplicated(out$datetimestamp))){
    
    out <- out[!duplicated(out$datetimestamp),]  
    
  }
  
  # remove rows with no datetimestamp
  out <- out[!is.na(out$datetimestamp), ]
  
  # convert output to data frame
  # retain only relevant columns
  out <- data.frame(
    datetimestamp = out$datetimestamp,
    out[, names(out) %in% nms], 
    row.names = seq(1, nrow(out))
    )
  
  # names as lower case
  names(out) <- tolower(names(out))

  # remove date from station_code, convert to swmpr class
  station_code <- gsub('[0-9]*$', '', station_code)
  out <- swmpr(out, station_code)
  
  if(trace) cat('\n\nData imported...')
  
  # return data frame
  return(out)
    
}

######
#' Obtain metadata for all stations
#' 
#' Obtain a \code{\link[base]{data.frame}} of metadata for all SWMP stations.
#' 
#' @export
#' 
#' @return A \code{data.frame} of SWMP metadata
#' 
#' @concept retrieve
#' 
#' @details This function retrieves data from the CDMO web services.  The computer making the request must have a registered IP address.  Visit the CDMO web services page for more information: \url{http://cdmo.baruch.sc.edu/webservices.cfm}. This is the CDMO equivalent of \code{exportStationCodesXML}.
#' 
#' @examples
#' \dontrun{
#' 
#' ## retrieve metadata for all sites
#' site_codes()
#' 
#' }
site_codes <- function(){
  
  # access CDMO web services
  serv <- "http://cdmo.baruch.sc.edu/webservices2/requests.cfc?wsdl"

  # get all station codes
  reply <- httr::GET(
    serv,
    query = list(method = 'exportStationCodesXMLNew'),
    )

  # parse reply from server
  out <- parser(reply)
  out$params_reported <- tolower(out$params_reported)
  
  # return output
  return(out)
  
}

######
#' Obtain metadata for a single reserve
#'
#' Get metadata for all the stations at a single SWMP reserve
#' 
#' @param nerr_site_id chr string of site, three letters
#' 
#' @concept retrieve
#' 
#' @export
#' 
#' @return An abbreviated \code{data.frame} of the SWMP metadata for the requested site
#' 
#' @details This function retrieves data from the CDMO web services.  The computer making the request must have a registered IP address.  Visit the CDMO web services page for more information: \url{http://cdmo.baruch.sc.edu/webservices.cfm}. This function is the CDMO equivalent of \code{NERRFilterStationCodesXMLNew}.
#' 
#' @examples
#' \dontrun{
#' 
#' ## retrieve metadata for all stations at a site
#' site_codes_ind('apa')
#' 
#' }
site_codes_ind <- function(nerr_site_id){
  
  # access CDMO web services
  serv <- "http://cdmo.baruch.sc.edu/webservices2/requests.cfc?wsdl"

  # get all station codes
  reply <- httr::GET(
    serv,
    query = list(  
      method = 'NERRFilterStationCodesXMLNew',
      NERRFilter = nerr_site_id
    )
  )

  # parse reply from server
  out <- parser(reply)
  out$params_reported <- tolower(out$params_reported)
  
  # return output
  return(out)
  
}
    