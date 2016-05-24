# Parse the specified input file and return the doc object
parseXML <- function(infile)
{
   xmlTreeParse(infile, trim=FALSE, addAttributeNamespaces=TRUE)
}

# Extract the top level node from an XML doc object
getTopNode <- function(doc)
{
   if (is.character(doc))
     stop('need to parse XML document with parseXML function')

   doc$doc$children[[1]]
}
