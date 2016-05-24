# Clean spaces and special characters from strings
clean_string <- function(x) {
  gsub(" |&|'|-|\\.", "", x)
}



# Regroup with characters
group_by_char <- function(x, vars) {
  dots <- vars %>%
    as.list %>%
    lapply(as.symbol)
  group_by_(x, .dots = dots)
}

#' Get list of valid columns
#'
#' List of valid columns accepted in \code{\link{query_master}}, \code{\link{sum_master}} and related functions.
#'
#' @seealso \code{\link{query_master}}, \code{\link{sum_master}}
#'
#' @export
valid_columns <- function() c("collection", "property", "name", "parent", "category", "region", "zone",
                              "period_type_id", "band", "sample", "timeslice", "time")


#' Test if elements in sample column are statistics
#'
#' In stochastic simulations, PLEXOS will return sample results and their statistics together. This function
#' makes it easy to separate them with a filter.
#'
#' @param x Vector of sample values from an rplexos query
#'
#' @examples
#' \dontrun{db <- plexos_open()}
#' \dontrun{res <- query_month(db, "Generator", "Generation")}
#' \dontrun{res %>% filter(sample_stats(sample))    # To obtain statistics}
#' \dontrun{res %>% filter(!sample_stats(sample))   # To obtain sample results}
#'
#' @export
is_sample_stats <- function(x)
  x %in% c("Max", "Min", "Mean", "StDev")

#' Get list of folders in the working directory
#'
#' List of existing folders in the working directory. This function is used when the wildcard symbol (\code{"*"})
#' is provided to the \code{\link{process_folder}} and \code{\link{plexos_open}} functions.
#'
#' @seealso \code{\link{setwd}}, \code{\link{process_folder}}, \code{\link{plexos_open}}
#'
#' @export
list_folders <- function() {
  f <- dir()
  f[file.info(f)$isdir]
}


#### Validation rules ####

# Check that object is valid rplexos databasae
check_rplexos <- function(x) {
  if(!inherits(x, "rplexos"))
    stop("db is not a valid database object. It should be created with plexos_open().", call. = FALSE)
}

# Delete file and give error if unsuccesfull
stop_ifnot_delete <- function(x) {
  # Error if file cannot be removed
  suppressWarnings(did.remove <- file.remove(x))
  if (!did.remove)
    stop("Unable to delete file: ", x, call. = FALSE)
}

# Check that a vector of characters are folder names
check_is_folder <- function(x) {
  if ((length(x) == 1L) && identical(x, "*")) {
    test <- TRUE
  } else {
    test <- all(file.exists(x)) & all(file.info(x)$isdir, na.rm = FALSE)
  }
  
  if (!test)
    stop(paste0("'folders' must be a vector of existing folders or the wildcard \"*\""), call. = FALSE)
}
