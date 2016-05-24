if(.Platform$OS.type == "windows") {
    .onLoad <- function(libname, pkgname) {
        ## The user might already have set this.
        co <- getOption("RCurlOptions")
        if (!"cainfo" %in% names(co)) {
            ## Might check env variable CURL_CA_BUNDLE, which is normally
            ## set on Windows for R >= 3.2.0
            co$cainfo <- file.path(libname, pkgname, "etc", "ca-bundle.crt")
            options(RCurlOptions = co)
        }
    }
}

