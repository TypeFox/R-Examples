get_pinterest <-
function(links, sleep.time=0) {
    pin.response <- data.frame()
    pin.call <- paste0("http://api.pinterest.com/v1/urls/count.json?callback=receiveCount&url=",links,"&format=json")
    if(.Platform$OS.type == "windows") { if(!file.exists("cacert.perm")) download.file(url="https://curl.haxx.se/ca/cacert.pem", destfile="cacert.perm") }
    if(.Platform$OS.type == "windows") { api_scrapper <- function(x) try(RCurl::getURL(x, cainfo = "cacert.perm", timeout = 240, ssl.verifypeer = FALSE)) } else { 
        api_scrapper <- function(x) try(RCurl::getURL(x, timeout = 240, ssl.verifypeer = FALSE)) }
    Sys.sleep(sleep.time)
    pin.response <- try(data.frame(jsonlite::fromJSON(gsub("receiveCount\\(|\\)", "", api_scrapper(pin.call)))))
    return(pin.response)
}
