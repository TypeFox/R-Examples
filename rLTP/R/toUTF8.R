
##' Convert encoding of Chinese string to UTF-8.
##' 
##' @title Convert encoding of Chinese string to UTF-8.
##' @param cnstring A Chinese string vector.
##' @return Converted vectors.
##' @author Jian Li <\email{rweibo@@sina.com}>

toUTF8 <- function(cnstring)
{
	cnstring <- .verifyChar(cnstring)
	strenc <- Encoding(cnstring)[which(Encoding(cnstring) != "unknown")][1]
	if (!identical(strenc, "UTF-8")) {
		if(isUTF8(cnstring, TRUE)) {
			OUT <- cnstring
			Encoding(OUT) <- "UTF-8"
			return(OUT)
		} else if(isGB2312(cnstring, TRUE)) {
			strenc <- "gb2312"
			#} else if(isBIG5(cnstring, TRUE)) {
			#	strenc <- "big5"
		} else if(isGBK(cnstring, TRUE)) {
			strenc <- "GBK"
		} else {
			if (is.na(strenc)) strenc <- "GBK"
		}
	}
	OUT <- iconv(cnstring, strenc, "UTF-8")
	return(OUT)
}

