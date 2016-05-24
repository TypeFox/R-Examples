#' Read and combine raw OpenSesame data.
#'
#' Read and combine multiple raw data files that were collected with OpenSesame
#' (Mathot, Schreij, & Theeuwes, 2012). The files will be merged into one
#' \link{data.frame}.
#'
#' OpenSesame generally produces an output \code{.csv} file for each participant
#' in the experiment. This is handy during data collection, but for the analysis
#' it is often useful to combine many such files into a single
#' \code{data.frame}. This is the single task of the \code{read_opensesame}
#' function, which loads all files from a given directory and attempts to
#' combine them into a \code{data.frame}.
#'
#' \code{read_opensesame} provides a wrapper around \link{read_bulk} to load the
#' raw data files. After loading, the different data files are merged using
#' \link[plyr]{rbind.fill}. This function can deal with varying column names
#' across files, and still places data into the appropriate columns. If a column
#' is not present in a specific file, it will be filled with \code{NA}.
#'
#' @inheritParams read_bulk
#'
#' @return A \link{data.frame} containing the merged raw data.
#'
#'   One column in the data.frame (\code{File}) contains the name of the raw
#'   data file. If the \code{subdirectories} option is set, an additional column
#'   (\code{Subdirectory}) with the name of the subdirectory is added.
#'
#' @references Mathot, S., Schreij, D., & Theeuwes, J. (2012). OpenSesame: An
#' open-source, graphical experiment builder for the social sciences.
#' \emph{Behavior Research Methods, 44}(2), 314-324.
#'
#' @seealso \link{read_bulk} for reading and combining multiple data files
#' that have other file formats.
#'
#' @examples
#' \dontrun{
#' # Read single raw data file from OpenSesame
#' raw_data <- utils::read.csv("raw_data/subject-1.csv",encoding = "UTF-8")
#'
#' # Merge all files in the main folder "raw_data"
#' # (which is in the current working directory)
#' raw_data <- read_opensesame(directory = "raw_data")
#'
#' # Merge files with file extension ".csv"
#' raw_data <- read_opensesame(directory = "raw_data",
#'   extension = ".csv")
#'
#' # Merge all files stored in separate folders
#' # within the folder "raw_data"
#' raw_data <- read_opensesame(directory = "raw_data",
#'   subdirectories = TRUE)
#'
#' # Merge all raw data stored in the folders "Session1"
#' # and "Session2" within the folder "raw_data"
#' raw_data <- read_opensesame(directory = "raw_data",
#'   subdirectories = c("Session1","Session2"))
#'
#' # Export merged data to a file using write.table
#' write.table(raw_data, file = "raw_data.csv",
#'   sep=",", row.names = FALSE)
#' }
#' @export
read_opensesame <- function(directory=".",
  subdirectories=FALSE,
  extension=NULL,
  data=NULL) {

  return(read_bulk(
    directory=directory,
    subdirectories=subdirectories,
    extension=extension,
    data=data,
    fun=utils::read.csv,
    sep=",", header=TRUE, dec=".",
    stringsAsFactors = FALSE,
    encoding="UTF-8"
  ))

}
