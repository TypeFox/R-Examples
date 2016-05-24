
find_default_fm <- 
function(){
    switch( .Platform$OS.type
          , windows = "explorer.exe"
          , unix    = switch( Sys.info()['sysname']
                            , Darwin = "open"
                            , "xdg-open"
                            )
          )
}

#' Open the file manager
#' 
#' @param dir Directory to open.
#' @param manager the file manager to use.
#' 
#' @export
file_manager <- 
function( dir     = getwd()                                         #< Directory to open
        , manager = getOption("file.manager", find_default_fm())    #< The file manager to use.
        ){
    "Open the file manager."
    dir <- normalizePath(dir)
    suppressWarnings(system2(manager, dir, wait=FALSE, invisible=FALSE))
    #! @return called for the side effect.
}

#' @export
#' @describeIn file_manager alias for \code{file_manager}
explorer <- file_manager

