

#' Concatenate a vector of strings
#' @name strConcat
#' 
#' @param ss vector of strings
#' @param sep a length-1 string used as separator, default to ""
#' @return a string
#' @examples
#' \dontrun{
#' strConcat(letters)
#' strConcat(letters, " ")
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @docType methods
strConcat = function(ss, sep = "") {
	stopifnot(length(sep) == 1)
	paste(ss, collapse = sep)
}


#' String Representation of a character vector
#' 
#' 
#' @param ss character. 
#' @param print_out logical. Whether to print out the string representation.
#' @param single_quote Logical, whether to use single quote for wrap strings. Default to TRUE, when set to FALSE, double quote is used. 
#' @param start_with_c Logical, whether the representation should start with "c(", when set to FALSE, "(" is used. Default to TRUE.
#' @return character.
#' @examples 
#' \dontrun{
#' strVectorRepr(letters[1:3]) == 'c("a", "b", "c")'
#' strVectorRepr(
#'   as.character(1:3)) == 'c("1", "2", "3")'
#' all(eval(parse(text = strVectorRepr(as.character(1:3)))) == 
#'   c("1", "2", "3"))
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
strVectorRepr = function(ss, print_out = FALSE, single_quote = TRUE, start_with_c = TRUE) {
	delim = if(single_quote) "'" else '"'
	start = if(start_with_c) "c(" else "("
	ss = strConcat(
			c(
					start,
					strConcat(
							paste(delim, ss, delim, sep = ""), 
							", "),
					")"
			)
	)
	if(print_out) {
		message(ss)
	}
	ss
}




#' String representation of a character vector for SQLite consumption
#' 
#' Transform a character vector (e.g. \code{c("a", "b")} into a string representation 
#' that can be used in a SQLite query (e.g. "('a', 'b')"). 
#' 
#' @param vec character.
#' @param print_out logical. Print out the string representation when set to TRUE.
#' @param single_quote logical. Whether to use single quote for each element. Use double quote if set to FALSE. Default to TRUE.
#' 
#' @author Kaiyin Zhong
#' @export
strVectorSQLRepr = function(vec, print_out = FALSE, single_quote = TRUE) {
	if(single_quote) {
		quoted_strings = strConcat(paste("'", vec, "'", sep = ""), ",")
	} else {
		quoted_strings = strConcat(paste('"', vec, '"', sep = ""), ",")
	}
	res = strConcat(c("(", quoted_strings, ")"))
	if(print_out) {
		message(res)
	}
	res
}

#' String representation of a numeric vector for SQLite consumption
#' 
#' Transform a numeric vector (e.g. \code{c(1, 2)} into a string representation 
#' that can be used in a SQLite query (e.g. "(1, 2)"). 
#' 
#' @param vec numeric.
#' @param print_out logical. Whether to print out the string representation.
#' 
#' @author Kaiyin Zhong
#' @export
numVectorSQLRepr = function(vec, print_out = FALSE) {
	res = strConcat(c(
					"(", 
					paste(as.character(vec), collapse = ", "), 
					")"
			))
	if(print_out) {
		message(res)
	}
	res
}



#' Retrive SNP positions from hg19 database
#' @param snps A vector of SNP names
#' @param rm_underscore Remove irregular chromosome names
#' @return data.frame
#' 
#' @author kaiyin
#' @export
snpPosSNP138 = function(snps, rm_underscore = TRUE) {
	cmdStart = "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -NA     -e \"select chrom, chromStart, chromEnd, name from hg19.snp138 where name in "
	tmpfile = tempfile()
	cmdEnd = sprintf("\" hg19 > %s", tmpfile)
	snplist = strVectorRepr(snps, start_with_c = FALSE)
	cmd = sprintf("%s %s %s", cmdStart, snplist, cmdEnd)
	system(cmd)
	res = read.table(tmpfile, header = FALSE)
	colnames(res) = c("chrom", "chromStart", "chromEnd", "SNP")
	hasUnderscore = stringr::str_detect(res$chrom, "_")
	res[!hasUnderscore, ]
}


