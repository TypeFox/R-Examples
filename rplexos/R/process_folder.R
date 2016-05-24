#' Convert PLEXOS files to SQLite databases
#'
#' Functions to process PLEXOS solution and input files and dump them into SQLite databases
#' that are easier to read.
#'
#' \code{process_solution} converts a single zipped PLEXOS solution into a SQLite database.
#' The resulting database file has a simplified structure so data can be retrieved easily.
#' The temporary database contains the original data structure in the PLEXOS solution file.
#'
#' \code{process_input} converts a single PLEXOS input file (in XML form) solution into a
#' SQLite database. The database includes the raw tables and a series of views that better
#' organize the data.
#'
#' \code{process_folders} is used to process one or more folders. It searches each folder for
#' PLEXOS solution or input files and applies \code{process_solution} or
#' \code{process_input} to each one of them.
#'
#' The default folder is the working directory. If the wildcard \code{"*"} is provided, all
#' the folders in the working directory will be processed (the list of folders if provided by
#' the \code{\link{list_folders}} function).
#'
#' Do not rename the SQLite databases created with these functions. Other code expects
#' those filenames to remain unchanged.
#' 
#' @param folders Folder(s) to process (See details)
#' @param file Single PLEXOS solution or input file to process
#' @param keep.temp Should temporary databases be preserved?
#'
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Process the folder with the input file provided by rplexos
#' location2 <- location_input_rplexos()
#' process_folder(location2)
#' 
#' # Other examples
#' \dontrun{process_folder()}
#' \dontrun{process_solution("HiWind/Model WWSIS_c_RT_CoreB_M01_SC3 Solution.zip")}
#' \dontrun{process_input("WWSIS model.xml")}
#'
#' @export
process_folder <- function(folders = ".", keep.temp = FALSE) {
  # Check inputs
  stopifnot(is.character(folders), is.logical(keep.temp), length(keep.temp) == 1L)
  check_is_folder(folders)
  
  # Check for wildcard
  if (length(folders) == 1L) {
    if (identical(folders, "*")) {
      folders <- list_folders()
    }
  }
  
  # Check that folders exist
  are.dirs <- file.info(folders)$isdir
  are.dirs[is.na(are.dirs)] <- FALSE
  if(!all(are.dirs)) {
    not.dirs <- folders[!are.dirs]
    stop(paste(not.dirs, collapse = ", "), " are not valid folders", call. = FALSE)
  }
  
  # Function to list PLEXOS files in each folder
  plexos_list_files <- function(df) {
    filename <- list.files(df$folder, ".xml$|.XML$|.zip$|.ZIP$", full.names = TRUE)
    
    data_frame(type = ifelse(grepl(".xml$|.XML$", filename), "I", "S"),
               filename)
  }
  
  # Get database file names
  df <- data.frame(folder = folders,
                   stringsAsFactors = FALSE) %>%
    group_by(folder) %>%
    do(plexos_list_files(.))
  
  # Error if all folders were empty
  if (nrow(df) == 0L)
    stop("No input/solution files were found", call. = FALSE)
  
  # Check for folders without databases
  folder.missing <- setdiff(folders, df$folder)
  if (length(folder.missing) > 0L) {
    warning("No databases found in folder",
            ifelse(length(folder.missing) == 1L, ": ", "s: "),
            paste(folder.missing, collapse = ", "),
            call. = FALSE)
  }
  
  # Create new id for identification on screen
  df2 <- df %>%
    group_by(type) %>%
    arrange(filename) %>%
    mutate(id = paste0(type, 1:n()))
  
  rplexos_message("Found files:")
  for (i in 1:nrow(df2))
    rplexos_message("\t", df2$id[i], ":\t", df2$filename[i], sep = "")
  rplexos_message("")
  
  # Process files
  rplexos_message("Processing files:")
  
  if (!is_parallel_rplexos()) {
    df2 %>%
      group_by(id) %>%
      do(if(.$type == "I") {
        process_input(.$filename)
        data.frame()
      } else {
        process_solution(.$filename, keep.temp)
        data.frame()
      })
  } else {
    foreach(i = df2$id, .packages = c("dplyr", "rplexos", "DBI", "RSQLite")) %dopar% {
      df3 <- df2 %>% filter(id == i)
      if (df3$type == "I") {
        process_input(df3$filename)
      } else {
        process_solution(df3$filename, keep.temp)
      }
    }
  }
  
  invisible(TRUE)
}
