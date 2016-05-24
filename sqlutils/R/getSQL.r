#' Returns the query as a string with the parameters set.
#' 
#' @param query the query name.
#' @param ... SQL parameters.
#' @return the SQL string with parameters replaced.
#' @export
getSQL <- function(query=NULL, ...) {
	sql <- getSQLRaw(query)
	parmvals <- unlist(list(...))
	parms <- getParameters(query)
	notset <- parms[!parms %in% names(parmvals)]
	doc <- sqldoc(query)
	if(length(notset) > 0) {
		params <- doc$params
		for(v in notset) {
			if(!is.null(params) & length(params[params$param == v, 'default']) > 0 &
					!is.na(params[params$param == v, 'default'])) {
				val <- params[params$param == v, 'default']
				val <- eval(parse(text=val))
				parmvals = c(parmvals, val)
				names(parmvals)[length(parmvals)] <- v
				warning(paste("The ", v, ' parameter has not been set. Using the default value of ',
							  val, sep=''))
			} else {
				stop(paste("The ", v, 
					" parameter has not been set and no default value exists", sep=''))
			}
		}
	}
	if(length(parmvals)>0) {
		for(i in 1:length(parmvals)) {
			sql <- gsub(paste(":", names(parmvals)[i], ":", sep=''), 
						parmvals[i], sql)
		}
	}
	return(sql)
}
