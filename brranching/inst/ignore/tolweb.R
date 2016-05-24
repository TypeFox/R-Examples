get_tol <- function(searchterm) {
  baseurl <- "http://tolweb.org/onlinecontributors/app?service=external&page=xml/GroupSearchService&group="
  url <- paste(baseurl, searchterm, sep = "")
  tt <- getURL(url)
  ttp <- xmlRoot(xmlTreeParse(tt))
  ID <- as.character(xmlAttrs(ttp[[1]], name="ID"))
  return(ID)
}

library(httr)
tt <- GET('http://tolweb.org/onlinecontributors/app?service=external&page=xml/GroupSearchService&group=Bembidion')
stop_for_status(tt)
content(tt)

out <- GET('http://tolweb.org/onlinecontributors/app?service=external&page=xml/TreeStructureService&node_id=146766')
stop_for_status(out)
txt <- content(out, "text")
cat(txt)

library("xml2")
txt2 <- content(out, "text")
xml <- xml2::read_xml(txt2)
one <- xml2::xml_children(xml)[[1]]
