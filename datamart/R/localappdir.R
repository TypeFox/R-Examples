#' Platform-independent local app folder
#'
#' determines local app folder based on Sys.info()["sysname"]. 
#'
#' @param appname    subdir in local app folder, default "R"
#' @param create     logical, default TRUE, create folder if non-existent
#'
#' @return character 
#' @export
localappdir <- function(appname="R", create=TRUE) {
    # resolve settings folder
    folder <- switch(tolower(Sys.info()["sysname"]),
      unix=file.path(Sys.getenv('HOME'), paste(".", tolower(appname), sep="")),
      windows=file.path(
        Sys.getenv('LOCALAPPDATA', Sys.getenv('APPDATA', getwd())), 
        paste(toupper(substr(appname,1,1)), tolower(substring(appname,2)), sep="")
      ),
      file.path(getwd(), appname)
    )
    
    # create settings folder if necessary
    if (create && !file.exists(folder)) dir.create(folder)

    return(folder)
}


