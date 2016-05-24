#' Import your individual ERP data files.
#'
#' \code{load.data} imports your individual ERP data files. File extensions must be .txt and
#' file names must be in the format: YourFile_Condition.txt (e.g., SS34_Positive.txt). Raw data files to be
#' imported should be organized as follows:
#' \itemize{
#'   \item each electrode must be a separate column
#'   \item voltages at each time point should be listed under the appropriate electrode
#'     column as rows
#'  \item no other data should be present in the raw data file (e.g., subject, condition,
#'    time, etc.)
#' }
#'
#' @param path The folder path containing your ERP files
#' @param condition In quotes, a string indicating which trial
#' type you will be importing (i.e., the condition indicated in the file name)
#' @param num.subs The number of files (subjects) to import for a given condition
#' @param epoch.st The earliest time point sampled in the ERP files, including
#'   the basline (e.g., -200)
#' @param epoch.end The final time point sampled in the ERP files
#' @param header Only accepts values of TRUE or FALSE. Used to specify whether or not there
#'   is an existing header row in the ERP files.  If there is no header, \code{load.data}
#'   will supply one (see details below).
#'
#' @details \itemize{
#'   \item Name each individual file following the format mentioned above (e.g., SS34_Positive.txt).
#'   \code{load.data} will ignore all text preceding the "_", and treat all text following the "_"
#'   as the \code{condition}, (e.g., Positive).  Use only one "_" in the file name (i.e., to separate
#'   your own naming convention from the \code{condition}); using multiple "_" characters will lead to
#'   faulty importing of data.  The erp.easy convention for subjects is a capital "S" followed by the
#'   number corresponding to the order in which the file was loaded (e.g., S1, S2, etc.). Subjects will
#'   be loaded into the "Subject" column of the returned data frame.
#'
#'   \item If no header is present in the ERP files, one will be supplied, using the standard R
#'   convention of a capital "V" followed by increasing integers (e.g., V1, V2, V3). Use these
#'   automatically assigned column name values to refer to the electrodes (unless a header is provided
#'   in the raw data file).
#'
#'   \item Enter the starting time of the baseline, if present in your individual files, in
#'   \code{epoch.st} (e.g., -200).
#'
#'   \item Once the desired data frames have been loaded, they can be
#'   \href{http://www.statmethods.net/input/exportingdata.html}{exported} as a number of
#'   different file types.
#'
#'   \item The sample rate will be calculated for you, based on the starting (\code{epoch.st})
#'   and ending (\code{epoch.end}) time points of the recording epoch and the number of time
#'   points in a given condition (the number of rows in your file for each condition).
#'}
#'
#' @note While importing data must be done using a separate function call for each condition,
#'   it can be convenient to use R's native \code{rbind.data.frame()} command to bind
#'   several loaded conditions (variables) into a single data frame consisting of multiple
#'   conditions. All erp.easy functions will act on all conditions included in the data frame
#'   passed to the function. For example, if you'd like to see all conditions plotted, simply
#'   use \code{rbind.data.frame()} to make a single data frame to pass to an erp.easy plotting
#'   function, and you will see all added conditions plotted simultaneously in the same figure
#'   (as opposed to making separate data frames for each condition, then passing each data
#'   frame separately to a function).
#'
#' @return A single, concatenated data frame of all electrode data for all
#'   subjects organized into columns, with three added columns:
#'
#' \enumerate{
#'   \item "Subject" containing repeating subject names
#'   \item "Stimulus" containing repeating condition names (e.g., Neutral)
#'   \item "Time" containing a repeating list of timepoints sampled
#' }
#'
#' @author Travis Moore
#'
#' @examples
#' \dontrun{
#' # Importing data for a condition named "Neutral" (file names: "Sub1_Neutral.txt",
#' "Sub2_Neutral.txt", etc.)
#' neutral <- load.data(path = "/Users/Username/Folder/", condition = "Neutral",
#' num.subs = 20, epoch.st = -200, epoch.end = 899, header = FALSE)
#'
#' # Adding imported data named "positive" to the imported "neutral" data
#' combo <- rbind.data.frame(neutral, positive)
#' }


load.data <- function(path, condition, num.subs, epoch.st, epoch.end, header = FALSE) {

  # restores original working directory upon exit
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd(path)

  # loads in all .txt files from specified directory
  files = list.files(path, pattern = "*.txt")

  # searches for specified condition
  cond.files = grepl(paste("^[^_]+_", condition, sep = ""), files)

  # count number of TRUE values in logical vector
  trues <- sum(cond.files, na.rm = TRUE)

  if (trues < 1) {
    stop(paste("NO FILES FOUND FOR CONDITION ", condition, "!", sep = ""))
  }

  # gets files for specified condition only
  sorted.files = files[cond.files]

  # check for number of existing files and entered subjects is the same
  if (trues != num.subs) {
    stop(paste("NUMBER OF SUBJECTS ", "(", num.subs, ")", " AND ACTUAL FILES ",
               "(", trues, ")", "DIFFERS!"))
  }

  # reads files into R
  tables <- lapply(paste(path, sorted.files, sep = ""), read.table,
                 header = header)

  # bind list into data frame
  data.df = plyr::ldply(tables)
  sublist = vector("list")
  for (i in 1:num.subs) {
    sublist[[i]] = c(rep(paste("S", i, sep = ""), (nrow(data.df)/num.subs)))
  }
  sublist = data.frame(matrix(unlist(sublist), ncol = 1))
  all.times = seq(epoch.st, epoch.end, 1)
  number = round(length(all.times)/(nrow(data.df)/num.subs), digits = 0)
  sampled.times = seq(epoch.st, epoch.end, number)
  stimlist = c(rep(condition, nrow(data.df)))
  data.df1 = cbind.data.frame(sublist, stimlist, sampled.times, data.df)
  colnames(data.df1)[1:3] <- c("Subject", "Stimulus", "Time")

  setwd(oldwd)

  return(data.df1)

  }


