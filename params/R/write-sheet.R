#' @rdname read_sheet
#' @export
#' @importFrom utils write.table
#'
write_sheet <- function(x, file, ext, ...){
	if(missing(ext))
		ext <- file_ext(file)

	dir.create(dirname(file), recursive = TRUE, showWarnings=FALSE)

	if(ext %in% c("tsv", "txt", "conf", "def", "mat")){
		write.table(x = x, file = file, sep = "\t", row.names = FALSE, quote = FALSE, ...)

	}else if(ext=="csv"){
		write.table(x = x, file = file, sep = ",", row.names = FALSE, quote = FALSE, ...)

	}else if(ext=="xlsx"){
		if (!requireNamespace('openxlsx', quietly = TRUE)) {
			stop("openxlsx needed for this function to work. Please install it.",
					 call. = FALSE)
		}
		openxlsx::write.xlsx(x, file = file, colNames = TRUE, ...)

	}else{
		stop("Sorry write_sheet does not recognize this file format: ", ext,
				 " please use tsv, csv or xlsx")
	}

}
