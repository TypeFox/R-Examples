URL_BASE <- "http://api.yummly.com/v1/api"
URL_GET <- paste(URL_BASE, "recipe", sep = "/")
URL_SEARCH <- paste(URL_BASE, "recipes", sep = "/")
URL_META <- paste(URL_BASE, "metadata", sep = "/")

auth_cache <- new.env()

func_cache <- new.env()

flavours <- c("sweet", "meaty", "sour", "bitter", "sweet", "piquant")

.onLoad <- function(libname, pkgname) {# nocov start
    options("yummlyr" = list(
        'log'                   = NULL
    ))
    # try loading YUMMLY_APP_ID and YUMMLY_APP_KEY from environment variables for testing
    app_id <- Sys.getenv("YUMMLY_APP_ID")
    app_key <- Sys.getenv("YUMMLY_APP_KEY")
    save_yummly_credentials(app_id, app_key)
}# nocov end
