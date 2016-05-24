##' @importFrom httpcache uncached
pollProgress <- function (progress_url, wait=1) {

    ## Configure polling interval. Will increase by rate (>1) until reaches max
    max.wait <- 30
    increase.by <- 1.2

    starttime <- Sys.time()
    timeout <- crunchTimeout()
    timer <- function (since, units="secs") {
        difftime(Sys.time(), since, units=units)
    }
    status <- uncached(as.numeric(crGET(progress_url)$progress))
    while (status < 100 && timer(starttime) < timeout) {
        Sys.sleep(wait)
        status <- uncached(as.numeric(crGET(progress_url)$progress))
        wait <- min(max.wait, wait * increase.by)
    }

    if (status != 100) {
        halt('Your process is still running on the server. It is currently ',
            status, '% complete. Check `uncached(crGET("',
            progress_url, '"))` until it reports 100% complete')
    }
    return(status)
}
