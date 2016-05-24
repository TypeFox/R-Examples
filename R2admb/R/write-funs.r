#'Write parameter and data files for ADMB
#'
#'Given base filenames and lists, write output files for starting parameter
#'values and data in a format suitable for input by AD Model Builder
#' from glmmADMB, by Hans Skaug
#'
#'@usage write_pin(name,L)
#' 
#'       write_dat(name, L, append=FALSE)
#'@aliases write_pin write_dat dat_write
#'@export write_pin write_dat dat_write
#'@param name (character) the base name of the file
#'@param L a list of objects to be written to file
#'@param append (logical) append to existing file?
#'@return Returns nothing; creates files in the current working directory as a
#'side effect
#'@note numeric vectors and matrices are the only objects that can be written
#'(at present)
#'@author Hans Skaug
#'@seealso \code{\link{read_pars}}
#'@keywords misc
write_pin <- "pin_write" <-
		function (name, L) 
{
	n <- nchar(name)
	if (substring(name, n - 3, n) == ".pin") 
		file_name <- name
	else file_name <- paste(name, ".pin", sep = "")
	cat("# \"", name, ".pin\" produced by pin_write() from R2admb ", 
			date(), "\n", file = file_name, sep = "")
	for (i in 1:length(L)) {
		x <- L[[i]]
		if (data.class(x) == "numeric") 
			cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
					append = TRUE)
		if (data.class(x) == "matrix") {
			cat("#", names(L)[i], "\n", file = file_name, append = TRUE)
			write.table(L[[i]], col.names = FALSE, row.names = FALSE, quote = FALSE, 
					file = file_name, append = TRUE)
			cat("\n", file = file_name, append = TRUE)
		}
	}
}


if (FALSE) {
	## test: can we read all ADMB examples without crashing?
	dir <- "/usr/local/src/admb/examples/admb/"
	dir <- "/usr/local/src/admb/examples/admb-re/"
	setwd(dir)
	## omit files with '.' (happen to be non-directories)
	L <- list.files(pattern="^[a-zA-Z_]+$")
	source("/home/ben/lib/R/pkgs/r2admb/pkg/R/admb-funs.R")
	for (i in seq_along(L)) {
		setwd(file.path(dir,L[i]))
		tpls <- gsub(".tpl","",list.files(pattern=".tpl"))
		for (j in seq_along(tpls)) {
			cat(L[i],tpls[j],"\n")
			invisible(read_tpl(tpls[j])$info)
		}
	}
}

## from glmmADMB, by Hans Skaug
write_dat <- "dat_write" <-
		function (name, L, append=FALSE) 
{
	n <- nchar(name)
	file_name <- if (tools::file_ext(name) == ".dat") {
				name
			} else paste(name, "dat", sep = ".")
	cat("# \"", file_name,"\" produced by dat_write() from R2admb ", 
			date(), "\n", file = file_name, sep = "", append=append)
	for (i in 1:length(L)) {
		x <- L[[i]]
		dc <- data.class(x)
		if (dc=="numeric") {
			cat("#", names(L)[i], "\n", L[[i]], "\n\n", file = file_name, 
					append = TRUE)
		} else {
			if (dc == "matrix") {
				cat("#", names(L)[i], "\n", file = file_name, append = TRUE)
				write.table(L[[i]], , col.names = FALSE, row.names = FALSE, quote = FALSE, 
						file = file_name, append = TRUE)
				cat("\n", file = file_name, append = TRUE)
			} else {
				stop(paste("can't handle data type '",dc,"' (variable ",names(L)[i],")",sep=""))
			}
		}
	}
}






