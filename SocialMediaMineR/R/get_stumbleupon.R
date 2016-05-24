get_stumbleupon <-
function(links, sleep.time=0) {
    stu.response <- data.frame()
    stu.call <- paste0("http://www.stumbleupon.com/services/1.01/badge.getinfo?url=",links)
    if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
    if(.Platform$OS.type == "windows") { api_scrapper <- function(x) try(RCurl::getURL(x, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE)) } else { 
        api_scrapper <- function(x) try(RCurl::getURL(x, timeout = 240, ssl.verifypeer = FALSE)) }
    Sys.sleep(sleep.time)
    stu.response <- try(data.frame(jsonlite::fromJSON(api_scrapper(stu.call))))
    return(stu.response)
}
