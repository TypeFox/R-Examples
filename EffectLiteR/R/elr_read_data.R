
#' Read Data File
#' 
#' Tries to determine the format of the data by the file ending and
#' chooses the appropriate function to read data. Currently supports .csv,
#' .dat, .txt, .sav, and .xpt and calls \code{\link[utils]{read.csv}},
#' \code{\link[utils]{read.table}}, \code{\link[foreign]{read.spss}},
#' \code{\link[foreign]{read.xport}} accordingly. The default values for
#' arguments depend on the function used to read data.
#' 
#' @param file Name of the file to read.
#' @param name Pure file name (without path to file) to read. 
#' If \code{file} includes a lengthy path name with many special characters, 
#' specifying this argument in addition to \code{file} may help the function 
#' to find the file ending.
#' @param header See \code{\link[utils]{read.table}}.
#' @param sep See \code{\link[utils]{read.table}}.
#' @param dec See \code{\link[utils]{read.table}}.
#' @param use.value.labels See \code{\link[foreign]{read.spss}}.
#' @param na.strings See \code{\link[foreign]{read.spss}}.
#' @return Object of class \code{"data.frame"}.
#' @export
elrReadData <- function(file, name=NULL, header="default", sep="default",
                        dec="default", use.value.labels="default",
                        na.strings="NA"){
  
  ptn <- "\\.[[:alnum:]]{1,5}$"
  if(is.null(name)){
    suf <- tolower(regmatches(file, regexpr(ptn, file)))
  }else{
    suf <- tolower(regmatches(name, regexpr(ptn, name))) 
  }
  
  ## convert arguments from shiny ui
  if(header=="yes"){header <- TRUE}
  if(header=="no"){header <- FALSE}
  if(sep=="semicolon"){sep <- ";"}
  if(sep=="white space"){sep <- ""}
  if(dec=="decimal point"){dec <- "."}
  if(dec=="decimal comma"){dec <- ","}
  if(use.value.labels=="yes"){use.value.labels <- TRUE}
  if(use.value.labels=="no"){use.value.labels <- FALSE}
  
  if(suf == ".csv"){
    if(header=="default"){header <- TRUE}
    if(sep=="default"){sep <- ","}
    if(dec=="default"){dec <- "."}
    return(read.csv(file, header=header, sep=sep, dec=dec,
                    na.strings=na.strings))
    
  }else if(suf == ".txt" || suf==".dat"){
    if(header=="default"){header <- FALSE}
    if(sep=="default"){sep <- ""}
    if(dec=="default"){dec <- "."}
    return(read.table(file, header=header, sep=sep, dec=dec,
                      na.strings=na.strings))
    
  }else if(suf == ".sav"){
    if(use.value.labels=="default"){use.value.labels <- TRUE}
    return(foreign::read.spss(file, to.data.frame=TRUE,
                              use.value.labels=use.value.labels))
    
  }else if(suf == ".xpt"){
    return(foreign::read.xport(file))
  }
  
}

