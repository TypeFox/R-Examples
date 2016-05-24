# Rd utility functions
# 
# Author: Renaud Gaujoux
# Created: Mar 25, 2013
###############################################################################

#getRdFile <- function(topic, package=NULL){
#	help_call <- substitute(help(topic, package = package, try.all.packages = TRUE), 
#			list(topic = topic, package = package))
#	
#	eval(help_call)
#}

# Borrowed from tools:::RdTags
RdTags <- function (Rd) 
{
    res <- sapply(Rd, attr, "Rd_tag")
    if (!length(res)) 
        res <- character()
    res
}


#' @importFrom tools Rd_db 
getRdTag <- function(topic, tag, package){
	# extract topic
	#rd <- utils:::.getHelpFile(file=getRdFile(topic, package=package))
    rd <- Rd_db(package=package)
    found <- FALSE
    i <- sapply(rd, function(x){
        if( found ) return()
        tags <- RdTags(x)
        w <- which(tags == "\\alias")
        if( length(w <- which(sapply(x[w], function(a) a[[1]] == topic))) ){
            found <<- TRUE
            rd
        }else NULL
    })
    if( !found ) stop("Could not find topic '", topic, "' in package '", package, "'")
    w <- which(!sapply(i, is.null))
    topic_rd <- rd[[w[1L]]]
    tags <- RdTags(topic_rd)
    if( !length(w <- which(tags == tag)) )
        stop("Could not find tag '", tag, "' in help topic ", package, "::", topic)
    topic_rd[w]
}

#' Format Rd Sections into LatTeX
#' 
#' This function extract sections from Rd files and convert them into 
#' LaTeX code. 
#' This can be useful to include Rd text into vignettes, hence keeping them 
#' up to date.
#' 
#' @section Example section: This is a nice section, with a bullet list: 
#' \itemize{
#' \item tata
#' \item toto
#' }
#' 
#' @param topic Rd topic
#' @param package package in which to search the topic
#' @param i index of the section to format
#' @param notitle logical that indicates if the section's title should be removed
#' 
#' @export
#' @examples
#' RdSection2latex('RdSection2latex', package = 'pkgmaker')
#' 
RdSection2latex <- function(topic, package, i=1L, notitle=TRUE){
	rdsec <- getRdTag(topic, tag="\\section", package = package)
	if( !length(rdsec) ) return()
	ltx <- capture.output(tools::Rd2latex(rdsec[i], fragment=TRUE))
	if( notitle ){
		parts <- stringr::str_match(ltx, "\\{Section\\}")
		w <- which(!is.na(parts[, 1]))
		ltx <- ltx[seq(w[1]+1, tail(w, 1)-1)]
	}
	ltx <- paste(ltx, collapse="\n")
	# remove link commands
	ltx <- gsub("\\\\LinkA\\{([^}]+)\\}\\{([^}]+)\\}", "\\2", ltx)
	
	cat(ltx)
	invisible()
}

