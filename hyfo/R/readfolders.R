
readData_folder <- function(folderName, keyword) {
  
  folderName <- paste(folderName, keyword, sep = '\\')
  
  fileNames <- list.files(folderName, pattern = '*.csv', full.names = TRUE)
  
  if (length(fileNames)==0) {
    fileNames <- list.files(folderName, pattern = '.TXT', full.names = TRUE)
    if (length(fileNames)==0) stop('Wrong keyword, initial has to be Upper-case')
    
    data <- collectData_txt_anarbe(folderName, rangeWord = c('D?a       ', -1, 'M?x.      ', -5))
    rownames(data) <- NULL
  } else {
    data <- collectData_csv_anarbe(folderName)
  }
  
  return(data)
}


# @importFrom utils choose.dir
#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 

readData <- function(keyword, folderName) {
  message('This function is only windows based, if you are using windows platform (real
          operational system, not virtual machine), and want to use this function, please
          contact the author (xuyuanchao37@gmail.com) for the windows version.')
#   message('Choose the main folder that, in it, there are different folders representing different gauging stations,
#            all the gauging stations have precipitation data, some of them also have discharge data,
#            this function is to open different gauging folders and read the data, arragen them together.')
#   message('\n\n
#            new file is a list based file and needs to be read by dget()')
  
#   fileNames <- list.files(folderName, full.names = TRUE)
#   data <- lapply(fileNames, FUN = readData_folder, keyword = keyword)
#   
#   names <- sapply(c(1:length(data)), function(x) colnames(data[[x]])[2])
#   names(data) <- names
#   
#   fileName <- file.choose(new = TRUE)
#   dput(data, file = fileName)
#   
#   return(data)
}
