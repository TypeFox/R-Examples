#' Reading excel files from Tesouro Direto
#'
#' Reads files downloaded by \code{\link{download.TD.data}}
#'
#' @param dl.folder Folder with excel files from tesouro direto
#' @param maturity A character vector with the desired maturities (in format
#'   ddmmyy, e.g. "010116"). If set to NULL, will download all maturities of
#'   all asset codes
#' @param cols.to.import Columns (numeric) to import from excel files (open and
#'   check the columns of the excel file from tesouro direto for details)
#' @param col.names Names of columns in final data.frame (same size as
#'   cols.to.import)
#'
#' @return A dataframe with the data
#' @export
#'
#' @examples
#' # Downloads data from tesouro direto (only 1 file for simplicity)
#'
#' dl.folder ='TD Files'
#' download.TD.data(asset.codes = 'LTN_2015', dl.folder = dl.folder)
#'
#' my.df <- read.TD.files(dl.folder = dl.folder, maturity ='010116')
#'
read.TD.files <- function(dl.folder = 'TD Files',
                              maturity = NULL,
                              cols.to.import = c(1,2,4),
                              col.names =  c('ref.date','yield.bid','price.bid')){

  #require(xlsx)

  # Error checking

  if (!is.character(dl.folder)) stop('Argument dl.folder should be a character')
  if (!is.null(maturity)) if(!is.character(maturity)) stop('Argument maturity should be a character or NULL')
  if (!is.numeric(cols.to.import)) stop('Argument cols.to.import should be numeric ')

  if(!any(cols.to.import==1)){
    stop('The first column represents the dates in the excel files. You should include it in cols.to.import')
  }

  if(length(cols.to.import)!=length(col.names)){
    stop('The lenghts of cols.cols.to.import and col.names should be equal')
  }

  if(any(cols.to.import>6)){
    stop('The maximum possible value in cols.to.import is 6. Thats the number of columns in the excel files.')
  }

  if (!dir.exists(dl.folder)){
    stop(paste('Cant find folder', dl.folder))
  }
  cat('\nReading xls data and saving to data.frame')

  my.f <- list.files(dl.folder, full.names = T)

  if (length(my.f)==0){
    stop(paste('Cant file xlsx files at',dl.folder))
  }


  my.df <- data.frame()
  for (i.f in my.f){

    cat(paste('\n Reading File = ', i.f, sep = ''))

   # Use capture.output so that no message "DEFINEDNAME" is shown
   # Details in: https://github.com/hadley/readxl/issues/82#issuecomment-166767220

   utils::capture.output( sheets <- readxl::excel_sheets(i.f) )

    # find maturities
    if (!is.null(maturity)){

      idx <- logical(length = length(sheets))
      for (i.matur in maturity){

        idx.temp <- stringr::str_detect(string = sheets, pattern = i.matur)
        idx <- idx|idx.temp

      }

      sheets <- sheets[idx]

    }

    if (length(sheets)==0) cat(paste('\n    Cant find maturities in file ',i.f ))


    for (i.sheet in sheets) {

      cat(paste('\n    Reading Sheet ', i.sheet, sep = ''))

      # Read it with readxl (use capture.output to avoid "DEFINEDNAME:"  issue)
      # see: https://github.com/hadley/readxl/issues/111
      utils::capture.output( temp.df <- readxl::read_excel(path = i.f,sheet = i.sheet) )

      # make sure it is a dataframe (readxl has a different format as output)
      temp.df <- as.data.frame(temp.df) 
      
      # control for columns

      if (max(cols.to.import)>ncol(temp.df)){
        stop(paste0('In file ',i.f, ' the number of columns is lower than ',max(cols.to.import),
                    '. Please supply values of cols.to.import that make sense.'))
      }

      # filter columns to import
      temp.df <- temp.df[,cols.to.import]
      n.col <- ncol(temp.df)

      # fix for different format of xls files

      colnames(temp.df) <- col.names
      temp.df[, c(2:n.col)] <- suppressWarnings(data.frame(lapply(X = temp.df[, c(2:n.col)], as.numeric)))

      temp.df <- temp.df[stats::complete.cases(temp.df), ]


      # fix for data in xls (for some files, dates comes in a numeric format and for others as strings)

      if (stringr::str_detect(as.character(temp.df$ref.date[1]),'/')){

        dateVec <- as.Date(as.character(temp.df[ , 1]), format = '%d/%m/%Y')

      } else {

        dateVec <- as.Date(as.numeric(temp.df$ref.date)-2,  origin="1900-01-01")
      }

      temp.df$ref.date <- dateVec
      temp.df$asset.code <- i.sheet

      my.df <- rbind(my.df, temp.df)
    }


  }
  
  my.df <- my.df[stats::complete.cases(my.df), ]

  return(my.df)
}

