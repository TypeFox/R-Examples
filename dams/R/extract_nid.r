#' Extract desired data on dams from pre-processed NID data
#'
#' @param sample_only logical flag indicating the desire to get only a sample 
#' of the NID data (which comes with this package) or the entire dataset 
#' @export
#' @examples
#' # sample NID data, 100 records only
#' dams_sample <- extract_nid()
#' 
#' # entire NID data, all the 74000+ records from bitbucket.org/rationshop
#' \dontrun{
#' dams_all <- extract_nid(sample_only = FALSE)
#' }
#'
extract_nid <- function(sample_only = TRUE) {
  
  #check inputs
  if (!is.logical(sample_only)) {
    stop("sample_only has to be either TRUE or FALSE!")
  }
  
  if (sample_only) {
    # get sample data
    nid_sample <- NULL
    data(nid_sample, envir = environment())
    
    return (nid_sample)
    
  } else {    
    # get complete data from bitbucket
    # code based on three tips - 
    # RCurl example on https
    # http://stackoverflow.com/questions/19890633/r-produces-unsupported-url-scheme-error-when-getting-data-from-https-sites
    # https://answers.atlassian.com/questions/122394/url-to-bitbucket-raw-file-without-commits

    nid_cleaned <- NULL
    
    nid_url <- "https://bitbucket.org/rationshop/packages/raw/master/nid_cleaned.txt"
    if(url.exists(nid_url, ssl.verifypeer = FALSE)) {
      message("downloading data from bitbucket. might take a few moments...")
      nid_data <- getURL(nid_url, ssl.verifypeer = FALSE)    
      nid_cleaned <- read.csv(text = nid_data, header = TRUE, quote = "", as.is = TRUE, sep = "\t")
    } else {
      stop("URL for the complete NID data does not exist!")
    }
    
    return (nid_cleaned)
  }  
}
