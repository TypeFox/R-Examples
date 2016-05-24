

#' @title
#' Read/Write sheets (automatically detect the file type and work accordingly)
#'
#' @description
#'
#' Read/Write sheets (automatically detect the file type and work accordingly)
#'
#' write_sheet requires version 0.3.1.
#'
#' \itemize{
#' \item \strong{tsv, txt, conf, def}: assumed to be tab-delimited
#' \item \strong{csv}: assumed to be comma delimited
#' \item \strong{xlsx}: microsoft excel, uses openxlsx to read the sheet.
#'   Also, it removes extra columns which often creep into excel files.
#' }
#'
#' @param x read: path to a file, to be read. write: a data.frame
#' @param ext determined using file extention. Specifying will override
#'
#' @param id_column all rows which have this column as blank are skipped. See details.
#' @param ... passed onto read.xlsx of openxlsx, read.table or read.csv2 depending on the file type.
#'
#' @param start_row supplied to read.xlsx
#' @param sheet supplied to read.xlsx, index or name of the sheet to be read from excel file. See \link[openxlsx]{read.xlsx}
#'
#' @param header first line is header? See \link{read.table}
#'
#' @param file write: output file name.
#' @param verbose verbosity level.
#'
#' @details
#'
#' \strong{Note: for excel sheets}:
#' \itemize{
#'   \item If \code{id_column} is missing, default if first column
#'   \item If \code{sheet} is missing, it automatically reads the first sheet
#' }
#'
#' \strong{Some important default values for tsv and csv files}:
#'
#' \code{stringsAsFactors = FALSE,comment.char = '#', strip.white=TRUE, blank.lines.skip=TRUE}
#'
#'
#' @importFrom tools file_ext
#' @importFrom utils read.table
#'
#' @examples
#'
#' ## read a excel sheet
#' sheet = read_sheet(system.file("extdata/example.xlsx", package = "params"))
#'
#' ## read a comma seperated sheet
#' csv = read_sheet(system.file("extdata/example.csv", package = "params"))
#'
#' ## read a tab seperate sheet
#' tsv = read_sheet(system.file("extdata/example.tsv", package = "params"))
#'
#'
#' # write sheets -------
#'
#' ## write a comma seperated sheet
#' write_sheet(sheet, "example.csv")
#'
#' ## write a tab seperated sheet
#' write_sheet(sheet, "example.tsv")
#'
#' ## write an excel seperated sheet
#' write_sheet(sheet, "example.xlsx")
#'
#' @export
read_sheet <- function(x, id_column, start_row = 1, sheet = 1, ext, header=TRUE, verbose = FALSE,  ...){
	if(missing(ext))
		ext <- file_ext(x)

	if(ext %in% c("tsv", "txt", "conf", "def")){
		mat <- utils::read.table(x, as.is=TRUE, sep="\t", header=header, stringsAsFactors = FALSE,
														 comment.char = '#', strip.white=TRUE, blank.lines.skip=TRUE, ...)

	}else if(ext=="csv"){
		mat <- utils::read.csv2(x, as.is=TRUE, sep=",", header=header, stringsAsFactors = FALSE,
														quote = "",
														comment.char = '#', strip.white=TRUE, blank.lines.skip=TRUE, ...)

	}else if(ext=="xlsx"){
		if (!requireNamespace('openxlsx', quietly = TRUE)) {
			stop("openxlsx needed for this function to work. Please install it.",
					 call. = FALSE)
		}
		mat <- openxlsx::read.xlsx(x, sheet = sheet, startRow = start_row, ...)
		message("Removing extra columns")
		mat = mat[,!grepl("^X", colnames(mat))]

	}else{
		stop("Sorry read_sheet does not recognize this file format: ", ext, " please use tsv, csv or xlsx (sheetname: sample_sheet)")
	}

	### ------ remove blank rows and columns
	if(missing(id_column)) {
		id_column = 1
		if(verbose)
			message("Reading file, using '", colnames(mat)[id_column], "' as id_column to remove empty rows.")
	}
	mat <- mat[!mat[, id_column] %in% c("", NA), ]
	return(mat)
}
