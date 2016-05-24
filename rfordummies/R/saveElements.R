
#' Saves a copy of the periodic table of elements as excel or csv file.
#' 
#' @param outfile File name
#' @param type Either excel or csv
#' 
#' @export
#' 
#' @examples 
#' saveElements(file.path(tempdir(), "elements.xlsx"))
#' saveElements(file.path(tempdir(), "elements.csv"), type = "csv")
#' list.files(tempdir(), pattern = "xlsx|csv", full.names = TRUE)
#' 
saveElements <- function(outfile, type = c("excel", "csv")){
  elements <- NA
  data("elements", package = "rfordummies", envir = parent.frame())
  type <- match.arg(type)
  switch(type, 
         excel = {
           if(!requireNamespace("XLConnect")) stop("Unable to create excel file. Install the XLConnect package and try again.")
           wb <- XLConnect::loadWorkbook(outfile, create = TRUE)
           XLConnect::createSheet(wb, name = "elements")
           XLConnect::writeWorksheet(wb, elements, sheet = "elements")
           XLConnect::saveWorkbook(wb)
         },
         csv = {
           write.csv(elements, file = outfile, row.names = FALSE)
         })
}

