##
## File:   ProcessXLSandXLSXFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the XLS and XLSX spreadsheet file processor.
##
processXLSandXLSXFile <- function(filePath = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(filePath, c("xls", "xlsx"))) return(NULL)

  # Read the workbook's worksheet names.
  worksheets <- readxl::excel_sheets(filePath)

  # Read the workbook's worksheets.
  workbook <- lapply(worksheets, readxl::read_excel, path = filePath)

  # Name the worksheets.
  names(workbook) <- worksheets

  # Return the results.
  workbook

}
