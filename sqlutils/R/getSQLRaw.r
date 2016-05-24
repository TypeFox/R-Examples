#' Returns the SQL from the file without the parameters replaced.
#' 
#' @param query the query name.
#' @return the unedited SQL statement.
getSQLRaw <- function(query) {
	f <- sqlFile(query)
    
	if(is.null(f)) { stop(paste("Cannot find query file for ", query, sep='')) }
	
	sql <- scan(f, what="character", sep=';', multi.line=FALSE, 
				comment.char=c("#"), quiet=TRUE, quote=NULL)
    
    sql <- ifelse( grepl("--", sql)
                  ,substr(sql, 0, regexpr("--", sql)-1 )
                  ,sql
                  )

	sql <- paste(sql, collapse=" ")

    sql <- gsub("  ", " ", sql)

	return(sql)
}
