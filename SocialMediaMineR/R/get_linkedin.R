get_linkedin <-
function(links, sleep.time=0) {
    lkn.response <- data.frame()
    lkn.call <- paste0("https://www.linkedin.com/countserv/count/share?url=",links,"&format=json")
    if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
    if(.Platform$OS.type == "windows") { api_scrapper <- function(x) try(RCurl::getURL(x, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE)) } else { 
        api_scrapper <- function(x) try(RCurl::getURL(x, timeout = 240, ssl.verifypeer = FALSE)) }
    Sys.sleep(sleep.time)
    lkn.response <- try(data.frame(jsonlite::fromJSON(api_scrapper(lkn.call))))
    return(lkn.response)
}
