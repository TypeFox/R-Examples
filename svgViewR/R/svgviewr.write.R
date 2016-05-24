svgviewr.write <- function(new_lines, save_as_fpath, append=FALSE){

	# IF FILE DOES NOT EXIST MAKE NEW FILE
	if(!file.exists(save_as_fpath)){write(new_lines, save_as_fpath);return(1)}

	# READ EXISTING FILE
	fileConn <- file(save_as_fpath)
	readLines.result <- readLines(fileConn)
	close(fileConn)

	if(append){
		# FIND LINE CONTAINING SVG_DOC TAG
		n_add <- grep('</svg_doc>', readLines.result)
		
		# COMBINE CURRENT AND NEW LINES
		new_lines <- c(readLines.result[1:(n_add-1)], new_lines, readLines.result[(n_add):length(readLines.result)])
	}else{
		# FIND LINES AT BEGINNING AND END OF SVG_DOC BLOCK
		n_start <- grep('<svg_doc xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" style=\"visibility:hidden;\" >', readLines.result)
		n_end <- grep('</svg_doc>', readLines.result)

		# COMBINE CURRENT AND NEW LINES
		new_lines <- c(readLines.result[1:(n_start)], new_lines, readLines.result[(n_end):length(readLines.result)])
	}

	# WRITE IN NEW LINES
	fileConn <- file(save_as_fpath)
	write(new_lines, fileConn)
	close(fileConn)
	
	1
}