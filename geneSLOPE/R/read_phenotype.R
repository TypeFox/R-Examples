#' Read phenotype from .fam file
#'
#' Reading phenotype data from file. It is assumed, that data is given in .fam file.
#' In this format, first column is family id (FID), second is individual id (IID),
#' third is Paternal individual ID (PAT), fourth is  Maternal individual ID (MAT),
#' fifth is SEX and sixth and last is PHENOTYPE.
#' If file has only four columns, then it is assumed that PAT and MAT columns are missing.
#' If there is only one column, then it is assumed that only phenotype is provided.
#'
#' @export
#' @param filename character, name of file with phenotype
#' @param sep character, field seperator in file
#' @param header logical, does first row of file contain variables names
#' @param stringAsFactors logical, should character vectors be converted to factors?
#'
#' @return object of class phenotypeData
#'
read_phenotype <- function(filename, sep=" ", header=FALSE, stringAsFactors=FALSE){
  if(!file.exists(filename))
    stop(paste("File", filename, "does not exist"))

  phe.data <- read.table(file = filename, header = header, sep = sep,
                         stringsAsFactors = stringAsFactors)
  if(ncol(phe.data)==6){
    y = phe.data[,6]
    y_info = phe.data
    colnames(y_info) = c("FID", "IID", "PAT", "MAT", "SEX", "PHE")
  } else if(ncol(phe.data)==4){
    y = phe.data[,4]
    y_info = phe.data
    colnames(y_info) = c("FID", "IID", "SEX", "PHE")
  } else if(ncol(phe.data)==1){
    y = phe.data[,1]
    y_info = NULL
  } else{
    stop("Incorrect data - see function description")
  }

  #returning phenotype data
  result <- structure( list(
    y = y,
    yInfo = y_info),
    class="phenotypeData")
  return(result)
}
