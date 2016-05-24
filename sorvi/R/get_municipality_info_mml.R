#' Get information of Finnish municipalities from Land Survey Finland.
#' (C) Maanmittauslaitos MML 2013. For details, see help(GetShapeMML).
#' @param ... Arguments to be passed
#' 
#' @return A data frame with municipality data
#' @export 
#' @references
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples \dontrun{tab <- get_municipality_info_mml()}
#' @keywords utilities

get_municipality_info_mml <- function (...) {

  # Load information table from Maanmittauslaitos
  map.id  <- "Yleiskartta-1000"
  data.id <- "HallintoAlue_DataFrame"

  # DO NOT SET sp <- NULL HERE; this will return NULL for the function
  # IN CONTRAST TO INTENDED OUTPUT!!!
  sp <- NULL	

  url <- paste(ropengov_storage_path(), "mml/rdata/", sep = "")
  filepath <- paste(url, map.id, "/", data.id, ".RData", sep = "")

  message(paste("Loading ", filepath, ". (C) MML 2013. Converted to RData shape object by Louhos. For more information, see https://github.com/avoindata/mml/", sep = ""))

  # Direct downloads from Github:
  # library(RCurl)
  # dat <- read.csv(text=getURL("link.to.github.raw.csv"))

  #load(url(filepath), envir = .GlobalEnv) # Returns a shape file sp
  load(url(filepath)) # Returns a data frame df

  # Vaasa and Hammarland have duplicated entries where only the 
  # enclave column differs. Remove that column, remove duplicated rows
  # and return the rest.
  df <- df[, -grep("Enklaavi", colnames(df))]
  df <- df[!duplicated(df), ]

  # Use harmonized municipality names as row.names
  # harmonized to match other data sets where 
  # slightly different versions of these names may be in use
  rownames(df) <- convert_municipality_names(df$Kunta.FI)

  # Order municipalities alphabetically
  df <- df[sort(rownames(df)), ]

  df

}
