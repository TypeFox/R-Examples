write.xml <- function(data, file=NULL) {
  if (!require(XML))
    stop("package XML must be installed")

  # Check that packages is available
#  if (! "XML" %in% row.names(installed.packages()) )
#    stop("package XML must be installed")

  if(is.null(file))
    stop("filename not specified")

  if (!is.data.frame(data))
    stop("data must be a data frame")

  # Start empty XML document tree
  doc <- XML::newXMLDoc()          
  # Start by adding a document tag at the root of the XML file
  root <- XML::newXMLNode("document", doc=doc)
  
  # Make output invisible
  invisible(
    # Iterate over all rows
    lapply(1:nrow(data),                 
           function(rowi) {
             r <- XML::newXMLNode("row", parent=root)   # Create row tag
             for(var in names(data)) {   # Iterate over variables
               XML::newXMLNode(var, data[rowi, var], parent = r)
             }
           }))            
  invisible(XML::saveXML(doc, file=file))
}
