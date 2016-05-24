#' Importing multiple file types and storing it into a Single R Dataframe
#'
#' importAllSdf function loads various data file types in a selected directory to dataframes,
#' combines all the data frames and stores it as a single data frame.
#' Note that all the files to be loaded should have the same number of columns for rbind to work
#'
#' @param path the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single dataframe containing all the files imported and rbinded into one.
#'
#'
#' @import  rChoiceDialogs
#'          tools
#'          readxl
#'          xlsx
#'          plyr
#'          utils
#'
#' @export

importAllSdf <- function(path=NULL)
{
  projectDirectory <- getwd();
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  if(files[[2]][[1]] == "xlsx") {
    importedDataFrame <- read_excel(toString(files[[1]][[1]]));
  }
  if(files[[2]][[1]] == "csv") {
    importedDataFrame <- read.csv(toString(files[[1]][[1]]),header = TRUE);
  }
  if(files[[2]][[1]]== "txt"){
    sample_text = read.delim2(toString(files[[1]][[1]]),header = TRUE)
    if(ncol(sample_text) == 1){
      sample_text = read.table(toString(files[[1]][[1]]), header = TRUE, sep = '|');
      if (ncol(sample_text) == 1) {
        sample_text = read.csv(toString(files[[1]][[1]]),header = TRUE, sep = '~');
        if (ncol(sample_text) == 1) {
          sample_text <- read.csv(toString(files[[1]][[1]]),header = TRUE, sep = '');
          importedDataFrame <- sample_text;
        }
        else {
          importedDataFrame <- sample_text;
        }
      }
      else {
        importedDataFrame <- sample_text;
      }
    }
    else {
      importedDataFrame <- sample_text;
    }
  }

  fileHeads <- names(importedDataFrame);

  for(i in 2:n){
    if(files[[2]][[i]] == "xlsx") {
      sample_text <- read_excel(toString(files[[1]][[i]]));
      names(sample_text) <- fileHeads;
      importedDataFrame <- rbind(importedDataFrame, sample_text);
    }
    else if(files[[2]][[i]] == "csv") {
      sample_text <- read.csv(toString(files[[1]][[i]]));
      names(sample_text) <- fileHeads;
      importedDataFrame <- rbind(importedDataFrame, sample_text);
    }
    if(files[[2]][[i]]== "txt"){
      sample_text <- read.delim(toString(files[[1]][[i]]), header = TRUE);
      if((ncol(sample_text) == 1)){
        sample_text <- read.table(toString(files[[1]][[i]]), header = TRUE, sep = '|');
        if ((ncol(sample_text) == 1)) {
          sample_text <- read.csv(toString(files[[1]][[i]]), header = TRUE, sep = '~');
          if ((ncol(sample_text) == 1)) {
            sample_text <- read.csv(toString(files[[1]][[i]]), header = TRUE, sep = '');
            names(sample_text) <- fileHeads;
            importedDataFrame <- rbind(importedDataFrame, sample_text);
            next;
          }
          else {
            names(sample_text) <- fileHeads;
            importedDataFrame <- rbind(importedDataFrame, sample_text);
            next;
          }
        }
        else {
          names(sample_text) <- fileHeads;
          importedDataFrame <- rbind(importedDataFrame, sample_text);
          next;
        }
        names(sample_text) <- fileHeads;
        importedDataFrame <- rbind(importedDataFrame, sample_text);
        next;
      }
      else {
        names(sample_text) <- fileHeads;
        importedDataFrame <- rbind(importedDataFrame, sample_text);
        next;
      }
    }
  }
  setwd(projectDirectory);
  return(importedDataFrame);
}

#' Importing multiple file types and storing it into a list of data frames
#'
#' importAllMdf function loads various data file types in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#' @export
importAllMdf <- function(path=NULL)
{
  projectDirectory <- getwd();
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]] == "xlsx") {
      sample_text <- read_excel(toString(files[[1]][[i]]));
      dflist[[toString(files[[1]][[i]])]] <- sample_text;
    }
    else if(files[[2]][[i]] == "csv") {
      sample_text <- read.csv(toString(files[[1]][[i]]));
      dflist[[toString(files[[1]][[i]])]] <- sample_text;
    }
    if(files[[2]][[i]]== "txt"){
      sample_text <- read.delim(toString(files[[1]][[i]]), header = TRUE);
      if((ncol(sample_text) == 1)){
        sample_text <- read.table(toString(files[[1]][[i]]), header = TRUE, sep = '|');
        if ((ncol(sample_text) == 1)) {
          sample_text <- read.csv(toString(files[[1]][[i]]), header = TRUE, sep = '~');
          if ((ncol(sample_text) == 1)) {
            sample_text <- read.csv(toString(files[[1]][[i]]), header = TRUE, sep = '');
            dflist[[toString(files[[1]][[i]])]] <- sample_text;
            next;
          }
          else {
            dflist[[toString(files[[1]][[i]])]] <- sample_text;
            next;
          }
        }
        else {
          dflist[[toString(files[[1]][[i]])]] <- sample_text;
          next;
        }
        dflist[[toString(files[[1]][[i]])]] <- sample_text;
        next;
      }
      else {
        dflist[[toString(files[[1]][[i]])]] <- sample_text;
        next;
      }
    }
  }

  setwd(projectDirectory);
  return(dflist);
}

#' Importing multiple csv files and storing it into a list of data frames
#'
#' importCsv function loads .csv file types in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#'
#' @export
importCsv <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]] == "csv") {
      sample_text <- read.csv(toString(files[[1]][[i]]));
      dflist[[toString(files[[1]][[i]])]] <- sample_text;
    }
  }
  setwd(projectDirectory);
  return(dflist);
}

#' Importing multiple csv files(";" delimited) and storing it into a list of data frames
#'
#' importCsv2 function loads .csv file types in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#'
#' @export
importCsv2 <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]] == "csv") {
      sample_text <- read.csv2(toString(files[[1]][[i]]));
      dflist[[toString(files[[1]][[i]])]] <- sample_text;
    }
  }
  setwd(projectDirectory);
  return(dflist);
}
#' Importing multiple excel files and storing it into a list of data frames
#'
#' importExcel function loads excel data file types(.xlsx, .xls) in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#' @param path  the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#'
#' @export
importExcel <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if((files[[2]][[i]] == "xlsx") || (files[[2]][[i]] == "xls")) {
      sample_text <- read_excel(toString(files[[1]][[i]]));
      dflist[[toString(files[[1]][[i]])]] <- sample_text;
    }
  }
  setwd(projectDirectory);
  return(dflist);
}

#' Importing multiple pipe delimitted files and storing it into a list of data frames
#'
#' importPipe function loads various text files which use a Pipe delimitter(|) in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path  the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#' @export
importPipe <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]]== "txt"){
      sample_text <- read.table(toString(files[[1]][[i]]), header = TRUE, sep = '|');
      if(ncol(sample_text)>1){
        dflist[[toString(files[[1]][[i]])]] <- sample_text;
      }
    }
  }
  setwd(projectDirectory);
  return(dflist);
}

#' Importing multiple tab delimitted files and storing it into a list of data frames
#'
#' importTab function loads various text files which uses a tab delimitter in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path  the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#'
#' @export
importTab <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]]== "txt"){
      sample_text <- read.delim(toString(files[[1]][[i]]), header = TRUE);
      if(ncol(sample_text)>1){
        dflist[[toString(files[[1]][[i]])]] <- sample_text;
      }
    }
  }
  setwd(projectDirectory);
  return(dflist);
}

#' Importing multiple tilde delimitted files and storing it into a list of data frames
#'
#' importTilde function loads various text files which uses a tilde delimitter(~) in a selected directory to separate dataframes,
#' and stores them as a list. Data frames can be accessed as list elements by using "listname$filename" or "listname[]".
#'
#'
#' @param path  the directory in which the files are stored.
#'        If path is not given, interactive dialog box will be used to select directory
#'
#' @return a single list of dataframes containing all the files imported and stored as dataframes.
#'
#' @import     rChoiceDialogs
#'             tools
#'             readxl
#'             xlsx
#'             plyr
#' @export

importTilde <- function(path=NULL)
{
  projectDirectory <- getwd();
  # Set the working directory for importing files.
  if(is.null(path)) {
    path <- rchoose.dir(default = getwd(), caption="Select the directory with the files to import");
  }
  setwd(path);
  fileNames <- list.files();
  filesExtensions <- file_ext(fileNames);
  files <- data.frame(fileNames,filesExtensions);
  n <- nrow(files);
  dflist = list();
  for(i in 1:n){
    if(files[[2]][[i]]== "txt"){
      sample_text <- read.csv(toString(files[[1]][[i]]), header = TRUE, sep = '~');
      if(ncol(sample_text)>1){
        dflist[[toString(files[[1]][[i]])]] <- sample_text;
      }
    }
  }
  setwd(projectDirectory);
  return(dflist);
}
