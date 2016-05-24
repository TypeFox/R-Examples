readHTMLLinks = getHTMLLinks =
function(doc, externalOnly = TRUE, xpQuery = "//a/@href", baseURL = docName(doc),
          relative = FALSE)
{
  if(is.character(doc))
     doc = htmlParse(doc)

    # put a . in front of the xpQuery if we have a node rather than a document.
  if(is(doc, "XMLInternalNode") && grepl("^/", xpQuery))
     xpQuery = sprintf(".%s", xpQuery)

  links = as.character(getNodeSet(doc, xpQuery))
  links = if(externalOnly)
             grep("^#", links, value = TRUE, invert = TRUE)
          else
             links

        #XXX Put base URL onto these links, relative!
  if(relative)
    sapply(links, getRelativeURL, baseURL)
  else
    links
}



getHTMLExternalFiles =
function(doc, xpQuery = c("//img/@src", "//link/@href", "//script/@href", "//embed/@src"),
           baseURL = docName(doc), relative = FALSE, asNodes = FALSE, recursive = FALSE)
{
  if(is.character(doc))
     doc = htmlParse(doc)

  if(asNodes)
     xpQuery = gsub("/@[a-zA-Z-]$+", "", xpQuery)

  nodes = getNodeSet(doc, xpQuery)

  if(asNodes)
    return(nodes)

  nodes = as.character(nodes)
  
  ans = if(relative)
          getRelativeURL(nodes, baseURL)
        else
          nodes
   # recursive.
  ans
}
