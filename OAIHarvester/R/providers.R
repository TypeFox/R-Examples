oaih_providers <-
local({
    .providers <- NULL
    function() {
        if(is.null(.providers)) {
            providers <-
                readHTMLTable("http://www.openarchives.org/Register/BrowseSites")
            providers <- providers[[2L]][3L : 5L]
            names(providers) <- c("name", "baseurl", "identifier")
            .providers <<- providers
        }
        .providers
    }
})
