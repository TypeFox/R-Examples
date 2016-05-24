##                    references.list                    ##
##      This code is part of the timetree package        ##
## copyright F.-S. Krah 2015 (last update: 2015-04-15)   ##
reference.list <- function() {
  url <- "http://www.timetree.org/reference_list.php"
  parse <- readHTMLTable(url)
  ref <- na.omit(parse[[1]])
  return(ref)
}