alrWeb <-
function (page = c("webpage", "errata", "primer"), script)
{
    script.page <- "http://www.stat.umn.edu/alr/Links/scripts/"
    page <- match.arg(page)
    urls <- c(webpage = "http://www.stat.umn.edu/alr/", 
        errata = "http://www.stat.umn.edu/alr/Links/errata.pdf", 
        primer = "http://www.stat.umn.edu/alr/Links/Rprimer.pdf")
    url <- urls[page]
    if(!missing(script)) url <- paste(script.page, script, ".R", sep="")
    browseURL(url)
}
