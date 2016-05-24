progressInit <- function(mode, n) {
    if (existingProgress())
        progressClose()
    if (! showingProgress() || n == 0)
        return()
    progressMessage(mode)
    con <- txtProgressBar(max = n, style = 3)
    setTxtProgressBar(con, 0)
    assign("progress", con, envir = .gridSVGEnv)
    assign("progressMode", mode, envir = .gridSVGEnv)
}

progressStep <- function(mode, n = 1) {
    if (!existingProgress() || !showingProgress() || diffProgressMode(mode))
        return()
    con <- get("progress", envir = .gridSVGEnv)
    curval <- getTxtProgressBar(con)
    setTxtProgressBar(con, curval + n)
    assign("progress", con, envir = .gridSVGEnv)
}

progressClose <- function() {
    if (! existingProgress() || ! showingProgress())
        return()
    con <- get("progress", envir = .gridSVGEnv)
    close(con)
    assign("progress", NULL, envir = .gridSVGEnv)
    assign("progressMode", NULL, envir = .gridSVGEnv)
}

progressMessage <- function(mode) {
    msg <- switch(mode,
                  grob = "Drawing grobs",
                  defs = "Drawing reference definitions",
                  pch = "Drawing pch definitions")
    message(msg)
}

existingProgress <- function() {
    exists("progress", envir = .gridSVGEnv) &&
    ! is.null(get("progress", envir = .gridSVGEnv))
}

diffProgressMode <- function(mode) {
    mode != get("progressMode", envir = .gridSVGEnv)
}

showingProgress <- function() {
    get("showProgress", envir = .gridSVGEnv)
}
