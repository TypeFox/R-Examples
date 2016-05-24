#' Produce an XML document of a table
#' 
#' Produce a document containing all data.table or data.frame rows 
#'
#' @param data  The data to be converted
#' @param name  The toplevel node name
#'
#' @return xml An XML object
#' 
#' @keywords XML xml data.frame data.table
#' @family helper
#' @export
#' 
#' @examples
#' df<-data.frame(nper=c(12,24),pmt=c(-10,-10),pv=c(110,220))
#' xml<-convertToXML(df,name='examples')
#'
#' @details
#' Code was taken from \url{https://stat.ethz.ch/pipermail/r-help/2010-February/228025.html} 
#' and amended, basic tests applied
#' 
#' @export
#'
#'

convertToXML <- function(data, name = "doc") {
    
    suppressWarnings(xml <- XML::xmlTree(name))
    
    xml$addTag(as.character(substitute(data)), close = FALSE)
    
    if (nrow(data) > 0) {
        for (i in 1:nrow(data)) {
            xml$addTag("row", close = FALSE)
            for (j in names(data)) {
                if (inherits(data, "data.table")) {
                  xml$addTag(j, data[i, j, with = FALSE])
                }
                
                if (!inherits(data, "data.table")) {
                  xml$addTag(j, data[i, j])
                }
            }
            xml$closeTag()
        }
    }
    xml$closeTag()
    return(xml)
} 
