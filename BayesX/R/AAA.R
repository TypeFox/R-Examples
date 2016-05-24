.BAYESX_CACHE <- new.env(FALSE, parent=globalenv())

.onLoad <- function(lib, pkg) {
    assign(".akimaStatus", FALSE, envir=.BAYESX_CACHE)
}

.onAttach <- function(lib, pkg) {
    startmessage <- paste("Note: Function plotsurf depends on akima which has\n",
          "a restricted licence that explicitly forbids commercial use.\n",
          "akima is therefore disabled by default and may be enabled by\n",
          "akimaPermit(). Calling this function includes your agreement to\n",
          "akima`s licence restrictions.\n",sep=" ")
    packageStartupMessage(startmessage, appendLF = FALSE)
}

.onUnload <- function(libpath) {
    rm(.BAYESX_CACHE)
}
