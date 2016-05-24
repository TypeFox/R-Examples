##                       getSYNS                     ##
##      This code is part of the rusda package       ##
##     F.-S. Krah 2015 (last update: 2015-07-11)     ## 

getSYNS <- function(x,  process, taxa)
{
    syns <- x[grep("\\s\r\n\t\t\t\t", x)]
    if (length(syns) > 0){
      syns <- unique(gsub(":\\s\r\n\t\t\t\t", "", syns))}
    else syns <- "nodata"
}