#' Open all PLEXOS databases
#' 
#' The default for \code{folders} is the working directory. If the wildcard \code{"*"} is provided, all
#' the folders in the working directory will be processed (the list of folders if provided by
#' the \code{\link{list_folders}} function).
#' 
#' Do not rename the SQLite databases created with the \code{\link{process_folder}} family of functions.
#' The code expects those filenames to remain unchanged.
#' 
#' @param folders character. Folder(s) where the data is located (each folder represents a scenario)
#' @param names character. Scenario names
#'
#' @seealso \code{\link{query_master}} to perform standard queries of data
#' @seealso \code{\link{query_sql}} to perform custom queries
#' 
#' @importFrom utils packageVersion compareVersion
#' @export
plexos_open <- function(folders = ".", names = folders) {
  # Check inputs
  stopifnot(is.character(folders), is.character(names))
  check_is_folder(folders)
  
  # Check for wildcard
  if (length(folders) == 1L) {
    if (identical(folders, "*")) {
      folders <- list_folders()
      names <- folders
    }
  }
  
  # Check that folder and names have the same length
  stopifnot(length(folders) == length(names))
  
  # Change default scenario name to something better than '.'
  if (length(folders) == 1L) {
    if (identical(folders, ".") & identical(names, ".")) {
      names <- "default"
    }
  }
  
  # Function to list PLEXOS files in each folder
  plexos_list_files <- function(df) {
    filename <- list.files(df$folder %>% as.character,
                           pattern = "rplexos.db$", full.names = TRUE)
    
    if (length(filename) == 0L)
      return(data.frame())

    data.frame(folder = df$folder,
               scenario = df$scenario,
               filename,
               stringsAsFactors = FALSE)
  }
  
  # Get database file names
  df <- data.frame(folder = folders,
                   scenario = factor(names, levels = names),
                   stringsAsFactors = FALSE) %>%
        rowwise() %>%
        do(plexos_list_files(.))
  
  # Error if all folders were empty
  if (nrow(df) == 0L)
    stop("No databases found in the list of folders.\n",
         "Did you forget to use process_folder()?",
         call. = FALSE)
  
  # Check for folders without databases
  folder.missing <- setdiff(folders, df$folder)
  if (length(folder.missing) > 0L) {
    warning("No databases found in folder",
            ifelse(length(folder.missing) == 1L, ": ", "s: "),
            paste(folder.missing, collapse = ", "),
            call. = FALSE)
  }
  
  # Open databases
  out <- df %>%
    ungroup() %>%
    mutate(position = 1:n()) %>%
    group_by(scenario, position, filename) %>%
    do(tables = get_list_tables(.$filename),
       properties = get_table(.$filename, "property")) %>%
    ungroup()
  
  # Add rplexos class to object
  class(out) <- c("rplexos", class(out))
  
  # Check the version of rplexos
  conf <- query_config(out)
  this.vers <- packageVersion("rplexos") %>% as.character
  if (!"rplexos" %in% names(conf)) {
    # rplexos is not even an entry in the config table
    warning("File(s) processed with an old version of rplexos. ",
            "Rerun process_folder() to avoid problems.",
            call. = FALSE)
  } else {
    # Compare to installed version
    comp <- sapply(conf$rplexos, compareVersion, this.vers)
    
    if (any(comp > 0)) {
      warning("File(s) processed with a newer version of rplexos. ",
              "Update rplexos or rerun process_folder() to avoid problems.",
              call. = FALSE)
    } else if (any(comp < 0)) {
      warning("File(s) processed with an older version of rplexos. ",
              "Rerun process_folder() to avoid problems.",
              call. = FALSE)
    }
  }
  
  out
}

# Create custom visualization for rplexos objects
#' @export
print.rplexos <- function(x, ...) {
  out <- x %>%
    group_by(scenario, position, filename) %>%
    summarize(tables = nrow(tables[[1]]),
              properties = nrow(properties[[1]])) %>%
    ungroup() %>%
    print()
}
