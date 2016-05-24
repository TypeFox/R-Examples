carWeb <-
function (page = c("webpage", "errata", "taskviews"), script, data)
{
    data.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/data/"
    script.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/scripts/"
    page = match.arg(page)
    urls = c(webpage = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/",
        errata = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/errata.html",
        taskviews = "http://cran.r-project.org/web/views")
    url <- urls[page]
    if(!missing(data)) {
       dfile <- unlist(strsplit(data, ".", fixed=TRUE))
       if(length(dfile) > 1) dfile <- dfile[1:(length(dfile)-1)]
       dfile <- paste(c(dfile, "txt"), collapse="." )
       url <- paste(data.page, dfile, sep="")}
    if(!missing(script)) {
       sfile <- unlist(strsplit(script, ".", fixed=TRUE))
       if(length(sfile) > 1) sfile <- sfile[1:(length(sfile)-1)]
       sfile <- paste(c(sfile, "R"), collapse="." )
       url <- paste(script.page, sfile, sep="")}   
    browseURL(url)
}