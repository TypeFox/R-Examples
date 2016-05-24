#' Collect data from csv for Anarbe case.
#' 
#' Collect data from the gauging stations in spain, catchement Anarbe
#' 
#' @param folderName A string showing the path of the folder holding different csv files.
#' @param output A boolean showing whether the output is given, default is T.
#' @return The collected data from different csv files.
#' @examples
#' 
#' #use internal data as an example.
#' file <- system.file("extdata", "1999.csv", package = "hyfo")
#' folder <- strsplit(file, '1999')[[1]][1]
#' a <- collectData_csv_anarbe(folder)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references 
#' 
#' \itemize{
#' \item http://meteo.navarra.es/estaciones/mapadeestaciones.cfm
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' @source http://meteo.navarra.es/estaciones/mapadeestaciones.cfm
#' @export
#' @importFrom utils tail
collectData_csv_anarbe <- function(folderName, output = TRUE){
  
  fileNames <- list.files(folderName, pattern='*.csv', full.names = TRUE)
  data <- lapply(fileNames, readColumn_csv_anarbe)
  data <- do.call('rbind', data)
  data <- data[, 1:2]
  data[, 1] <- as.Date(data[, 1], format = '%d/%m/%Y')
  
  #newFileName <- file.choose(new = T)
  #write.table(data_new,file=newFileName,row.names = F, col.names = F,sep=',')
  a <- unlist(strsplit(folderName, '\\\\|/'))
  tarName <- tail(a, 2)[1]
  colnames(data) <- c('Date', tarName)
  
  if (output) return(data)
}


readColumn_csv_anarbe <- function(fileName){
  data <- read.csv(fileName, skip = 4)
  endIndex <- which(data == '', arr.ind = TRUE)[1]-1
  
  data <- data[1:endIndex, ]
  
  if (!is.null(levels(data[, 2]))) {
    data[, 2] <- as.numeric(levels((data[, 2])))[data[, 2]]
  }
  
  colnames(data) <- c('Date', 'target')
  message(fileName)
  
  return(data)
}



#' Collect data from different excel files
#' 
#' @param folderName A string showing the folder path.
#' @param keyword A string showing the extracted column, e.g., waterLevel, waterBalance.
#' @param output A boolean showing whether the output is given.
#' @return The collected data from different excel files.
#' @export
#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
# @importFrom utils write.table
collectData_excel_anarbe <- function(folderName, keyword = NULL, output = TRUE){
  
  message('In order to make "hyfo" easier to be installed, this part is commented,
          check original R file in your computer or go to 
          https://github.com/Yuanchao-Xu/hyfo/blob/master/R/collectData_excel.R
          for ideas.')
  
  
  #   newFileName <- file.choose(new = TRUE)
  #   message ('new file should be located a different location than the excel folder, 
  #          in order to avoid error.
  #          At least 2 excels should be in the folder\n')
  #   
  #   message ('this function only applies to strange spain dem operation record file, and this strange file changes
  #          its format in the middle of the record. For other applications, some tiny changes needs to be made.')
  #   if (is.null(keyword)) stop('key word is needed, e.g."waterLevel".')
  #   
  #   fileNames <- list.files(folderName, pattern = '*.xls', full.names = TRUE)
  #   data <- lapply(fileNames, FUN = readColumn_excel_anarbe, keyword = keyword)
  #   checkBind(data, 'rbind')
  #   data <- do.call('rbind', data)
  #   
  #   data_new <- data.frame(data)
  #   
  #   data_new <- data_new[order(data_new[, 1]), ]
  #   
  #   
  #   startDate <- data_new[1, 1]
  #   endDate <- data_new[length(data_new[, 1]), 1]
  #   
  #   Date <- as.factor(seq(startDate, endDate, by = 1))
  #   
  #   if (length(Date) != length(data_new[, 1])) stop('check if the excel files are continuous')
  #   
  #   colnames(data_new) <- c('Date', keyword)
  #   
  #   write.table(data_new, file = newFileName,
  #               row.names = FALSE, col.names = TRUE, sep = ',')
  #   if(output == TRUE) return(data_new)
}

# 
# @importFrom xlsx read.xlsx
# readTable_excel_anarbe <- function(fileName){
#   
#   index <- tail(strsplit(fileName, '\\.|\\ ')[[1]], 3)
#   raw_year <- index[1]
#   raw_mon <- index[2]
#   
#   raw <- read.xlsx(fileName, sheetName='A')
#   startRow <- which(raw == 'COTA', arr.ind = TRUE)[1]+4
#   startCol <- which(raw == 'COTA',arr.ind = TRUE)[2]-1
#   stopRow <- which(raw =='TOTAL',arr.ind = TRUE)[1]-1
#   stopCol1 <- startCol + 17
#   stopCol2 <- which(raw == 'SUPERFICIE', arr.ind = TRUE)[2]
#   data <- cbind(raw[startRow:stopRow,startCol:stopCol1], raw[startRow:stopRow,stopCol2])
#   
#   
#   yearIndex <- rep(raw_year, stopRow-startRow+1)
#   monIndex <- rep(raw_mon, stopRow-startRow+1)
#   
#   data <- cbind(yearIndex, monIndex, data)
#   return(data)
# }
# # 
# @importFrom utils tail
# readColumn_excel_anarbe <- function(fileName, keyword = NULL){
#   
#   index <- tail(strsplit(fileName, '\\.|\\ ')[[1]],3)
#   year <- as.numeric(index[1])
#   mon <- as.numeric(index[2])
#   
#   if (year == 99) {
#     year = year + 1900
#   } else year = year + 2000
#   
#   word = c('COTA', 'Cota\n(m)', 'TOTAL', '  TOTAL')
#   
#   if (keyword == 'waterLevel') {
#     searchWord <- c('COTA', 'Cota\n(m)')
#   } else if (keyword == 'discharge_ERE') {
#     searchWord <- c('AF.ERE-', 'Caudal\n(m??/s)')
#   } else if (keyword == 'waterBalance') {
#     searchWord <- c('INCREMENTO', 'al Canal Bajo', 'AFORO',
#                     'Variaci??n\nvolumen embalsado')
#   } else if (keyword == 'surfaceArea') {
#     searchWord <- c('SUPERFICIE', 'SUPERFICIE')
#   } else if (keyword == 'volume') {
#     searchWord <- c('EMBALSADO', 'Volumen\n(m????)')
#   }
#   
#   
#   if (year == 1999 | year < 2009 | (year == 2009 & mon < 5)) {
#     raw <- xlsx::read.xlsx(fileName, sheetName = 'A')
#     startIndex <- which(raw == word[1], arr.ind = TRUE)
#     endIndex <- which(raw == word[3], arr.ind = TRUE)
#     startRow <- startIndex[1]+4
#     endRow <- endIndex[1]-1
#     
#     dayCol <- endIndex[2]
#     day <- raw[startRow:endRow, dayCol]
#     
#     targetCol <- which(raw == searchWord[1], arr.ind = TRUE)[2]
#     
#     if (is.na(targetCol)) stop(sprintf('capture nothing in %s', fileName))
#     
#     if (keyword == 'waterBalance') {
#       targetStart <- targetCol
#       targetEnd <- which(raw == searchWord[3], arr.ind = TRUE)[2]
#       a <- raw[startRow:endRow, targetStart:targetEnd]
#       a <- sapply(a, function(x) as.numeric(levels(x)[x]))
#       
#       if (year == 1999 & mon == 4) {
#         
#         target <- data.frame(a[, 2] * 86.4, a[, 5] * 86.4, rep(NA, dim(a)[1]), a[, 6] * 86.4,
#                              a[, 4] * 86.4, a[, 11] * 86.4, a[, 3], a[, 7], rep(NA, dim(a)[1]), a[, 1])
#       } else {
#         target <- data.frame(a[, 2] * 86.4, a[, 5] * 86.4, a[, 6] * 86.4, a[, 7] * 86.4, 
#                              a[, 4] * 86.4, a[, 12] * 86.4, a[, 3], a[, 8], rep(NA, dim(a)[1]), a[, 1])
#       }   
#       
#     } else {
#       target <- raw[startRow:endRow, targetCol]
#       if (keyword == 'discharge_ERE') target <- as.numeric(levels(target))[target]/1000
#     }
#     
#   } else {
#     raw <- read.xlsx(fileName,sheetName = 'parte del embalse')
#     startIndex <- which(raw == word[2], arr.ind = TRUE)
#     endIndex <- which(raw == word[4], arr.ind = TRUE)
#     startRow <- startIndex[1]+1
#     endRow <- endIndex[1]-2
#     
#     dayCol <- endIndex[2]
#     day <- raw[startRow:endRow, dayCol]
#     targetCol <- which(raw == searchWord[2], arr.ind=TRUE)[2]
#     if (is.na(targetCol)) stop(sprintf('capture nothing in %s', fileName))
#     
#     if (keyword == 'waterBalance') {
#       targetStart <- targetCol
#       targetEnd <- which(raw == searchWord[4], arr.ind=TRUE)[2]
#       target <- raw[startRow:endRow, targetStart:targetEnd]
#       
#     } else {
#       target <- raw[startRow:endRow, targetCol]
#     }
#     
#   }
#   
#   
#   startDate <- as.Date(paste(year, mon, day[1], sep = '-'))
#   endDate <- as.Date(paste(year, mon, tail(day,1), sep = '-'))
#   
#   Date <- seq(startDate, endDate, 1)
#   output <- data.frame(Date, as.vector(target))
#   colnames(output) <- c('Date', seq(1, dim(output)[2] - 1))
#   message(fileName)  
#   return(output)
#   
# }
# 





#' collect data from different txt.
#' 
#' @param folderName A string showing the folder path.
#' @param output A boolean showing whether the result is given.
#' @param rangeWord A list containing the keyword and the shift. 
#' defaut is set to be used in spain gauging station.
#' @examples
#'   
#' #use internal data as an example.
#' 
#' \dontrun{
#' file <- system.file("extdata", "1999.csv", package = "hyfo")
#' folder <- strsplit(file, '1999')[[1]][1]
#' a <- collectData_txt_anarbe(folder)
#' }
#'
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references 
#' 
#' \itemize{
#' \item http://www4.gipuzkoa.net/oohh/web/esp/02.asp
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' 
#' @source http://www4.gipuzkoa.net/oohh/web/esp/02.asp
#' @return The collected data from different txt files.
#' @export
#' @importFrom utils tail
collectData_txt_anarbe <- function(folderName, output = TRUE, rangeWord = c('Ene       ', -1, 
                                                                            'Total     ', -6)){
  #All the code should be ASCII encode, so there should be no strange symbol.
  if (is.null(rangeWord)) {
    stop('rangeWord consists of 4 elements:
         1. start word which program can recognise.
         2. shift1, the shift needs to be made. E.g. start word is in line 7, and program
         should read file from line 9, then shift is 9-7 = 2.
         3. end word, as start word
         4. shift2, same as shift1, sometimes can be negative
         
         E.g. rangeWord=c(\"aaa\",2,\"bbb\",-2)
         if no rangeWord, just input c(NULL,NULL,NULL,NULL)')
    
  }
  
  
  fileNames <- list.files(folderName, pattern = '*.TXT', full.names = TRUE)
  
  data <- lapply(fileNames, FUN = readColumn_txt_anarbe, rangeWord = rangeWord)
  
  data <- do.call('rbind', data)
  
  a <- unlist(strsplit(folderName, '\\\\|/'))
  tarName <- tail(a, 2)[1]
  colnames(data) <- c('Date', tarName)
  
  #newFileName <- file.choose(new = T)
  message('new file should be located a different location than the excel folder,
          in order to avoid error.
          At least 2 excels should be in the folder')
  
  #write.table(data_new,file=newFileName,row.names = F, col.names = F,sep=',')
  
  
  if (output == TRUE) return(data)
  
}  



anarbe_txt <- function(dataset, x1, x2){
  
  data <- as.matrix(dataset[x1:x2, 2:13])
  startYear <- data[1, 6]
  
  data <- data[5:35, ]
  
  date <- which(data != '          ', arr.ind = TRUE)
  startDate <- date[1, ]
  
  endDate <- date[length(date[, 1]), ]
  
  startDate <- as.Date(paste(startYear, startDate[2], startDate[1], sep = '-'))
  endDate <- as.Date(paste(startYear, endDate[2], endDate[1], sep = '-'))
  
  Date <- as.factor(seq(startDate, endDate, 1))
  
  dim(data) <- c(length(data), 1)
  
  data <- as.numeric(data[which(data != '          '), ])
  
  if (length(data) != length(Date)) {
    stop('check original txt file. for missing value, the symbol is "--", check
         if this symbol is missing somewhere')
  }
  
  output <- data.frame(Date = Date, target = data)
  
  return(output)
  }


#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' @importFrom utils read.fwf
readColumn_txt_anarbe <- function(fileName, keyword = NULL, rangeWord = NULL){
  
  a <- read.fwf(fileName, widths = rep(10,13))#read file with fixed width
  
  startRow <- which(a == rangeWord[1], arr.ind = TRUE)[, 1]
  startRow <- startRow + as.numeric(rangeWord[2])
  
  endRow <- which(a == rangeWord[3], arr.ind = TRUE)[, 1]
  endRow <- endRow + as.numeric(rangeWord[4])
  
  data <- mapply(FUN = function(x1, x2) anarbe_txt(dataset = a, x1, x2), startRow, endRow)
  
  data_new <- data.frame(Data = unlist(data[1, ]), target = unlist(data[2, ]))
  message(fileName)
  return(data_new)
}
