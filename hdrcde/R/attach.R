.onAttach <- function(...)
{
    library(help=hdrcde)$info[[1]] -> version
    version <- version[pmatch("Version",version)]
    um <- strsplit(version," ")[[1]]
    version <- um[nchar(um)>0][2]
    packageStartupMessage(paste("hdrcde",version,"loaded\n"))
}
