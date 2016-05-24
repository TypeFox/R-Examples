# Package load hook that set up the default options
.onLoad <- function(lib, pkg){
    sn <- tolower(Sys.info()["sysname"])
    term <- tolower(Sys.getenv("TERM"))
    if(term %in% c("ANSI", "apple_terminal")){
        options(style = "ANSI")
    } else if(sn == "windows" || term == "") {
        # RStudio also ends up here
        options(style = "off")
    } else {
        # Presumably linux or unix
        options(style = "xterm-256color")
    }
    style.palette("dark")
}

#' Color termnal output
#'
#' Talk about the different options. style.mode, style.palette.
#' 
#' @param ... Sent to \code{\link{cat}}.
#' @param fg Foreground color i.e. color of the text. Can be any number in
#'   [0, 255] or a string such as \code{"grey"} or \code{"dark red"} (for the
#'   basic 16 colors).
#' @param bg Background color. Takes the same values as \code{fg}.
#' @param font A vector containg any combination of the following font styles:
#'   \code{"normal"}, \code{"bold"}, \code{"underline"},
#'   \code{"blink"} (renders as bold on most terminals), \code{"inverse"}.
#'   Note that these may not be rendered on all terminals.
#' @param mode Escape code mode.
#' @return Nothing, sends all output to \code{\link{cat}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @seealso style
#' @examples
#' style("Blue suede shoes\n", bg="blue")
#' style(fg="red")
#' cat("everything is red now!")
#' style(NULL)
#' cat("but not anymore!")
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
style <- function(..., fg=NA, bg=NA, font=NA, mode=c("xterm-256color", "ANSI", "off")){
    if(is.null(getOption("style")))
        options(style=match.arg(mode))
    if(missing(mode)){
        mode <- getOption("style")
    } else if(missing(...)){
        options(style = mode)
        if(mode != "off") options(style.if.on = mode)
    }

    named.colors <- switch(mode,
        `ANSI` = c("black", "red", "green", "yellow", "blue", "magenta", "cyan", "white"),
        `xterm-256color` = c("black", "dark red", "dark green", "dark yellow",
                             "dark blue", "dark magenta", "dark cyan", "grey",
                             "dark grey", "red", "green", "yellow", "blue",
                             "magenta", "cyan", "white"),
        character(0)
    )
    if(is.null(fg)) fg <- NA
    if(is.null(bg)) bg <- NA
    if(is.null(font)) font <- NA
    if(is.character(fg)) fg <- match(fg, named.colors)-1
    if(is.character(bg)) bg <- match(bg, named.colors)-1

    fonts <- c(normal=0, bold=1, underline=4, underscore=4, blink=5, inverse=7)
    if(is.character(font)) font <- fonts[font]

    if(mode == "xterm-256color"){
        escape.str <- c(
            if(!is.na(fg)) sprintf("\033[38;5;%im", fg) else NULL,
            if(!is.na(bg)) sprintf("\033[48;5;%im", bg) else NULL
        )
    } else if(mode == "ANSI"){
        escape.str <- c(
            if(!is.na(fg)) sprintf("\033[%im", 30 + xterm256.to.ANSI(fg)) else NULL,
            if(!is.na(bg)) sprintf("\033[%im", 40 + xterm256.to.ANSI(bg)) else NULL
        )
    } else {
        escape.str <- ""
    }
    if(!is.na(font)) escape.str <- c(escape.str, sprintf("\033[%im", font))

    out <- paste(escape.str, collapse="")
    if(missing(...)){
        # Make the settings stick
        options(style.mode = mode)
    } else {
        out <- sprintf("%s%s\033[0m", out, paste(...))
    }
    structure(out, class="xtermStyle")
}
#' @export
#' @rdname style
style.clear <- function(...) style("", ...)
#' @export
#' @rdname style
style.off <- function()
    options(style = "off")
#' @export
#' @rdname style
style.on <- function()
    options(style = getOption("style.if.on"))

#' Print a style object
#'
#' @param x \code{\link{style}} object.
#' @param ... Ignored, kept for S3 consistency.
#' @seealso style
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
print.xtermStyle <- function(x, ...){
    cat(x)
}

#' Define the palette for auto-styling
#'
#' @param x Palette name or definition.
#' @param ... Settings that override the defaults.
#' @return Nothing, sets the approprite \code{\link{options}} variables.
#' @examples
#' style.palette("light")
#' style.palette(list(fg=c(numeric=10, character=13)))
#' style.palette("light", fg=c(numeric=10, character=13))
#' @seealso style.auto, display.xterm.colors
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
style.palette <- function(x, ...){
    if(is.character(x)){
        x <- get.palette(x)
    } else if(missing(x)){
        x <- getOption("style.palette", list())
    } else if(!is.list(x)){
        stop("Invalid palette")
    }
    if(!missing(...)){
        l <- list(...)
        stopifnot(all(sapply(l, is.vector)))
        x[names(l)] <- Map(function(x, y) c(x, y[!names(y) %in% names(x)]),
                           l, x[names(l)])
    }
    options(style.palette = x)
}
#' @rdname style.palette
#' @export
get.palette <- function(x = c("dark", "light")){
    switch(match.arg(x),
        dark = list(
            fg = c(
                null = 246,    # grey
                character = 2, # green
                numeric = 6,   # cyan
                factor = 5,    # magenta
                logical = 3,   # yellow
                list = 1,      # red
                `function` = 4 # blue
            ),
            bg = c(
                scalar = NA,
                vector = 255,
                matrix = 195,
                array = 224
            ),
            levels = xterm.pal()$paired[c(FALSE, TRUE)],
            range = xterm.pal()$DownUp,
            logical = xterm.pal()$Accent[c(5,3)]
        ),
        light = list(
            fg = c(
                null = 236,     # Dark grey
                character = 10, # green
                numeric = 14,   # cyan
                factor = 13,    # magenta
                logical = 11,   # yellow
                list = 9,       # red
                `function` = 33 # blue
            ),
            bg = c(
                scalar = NA,
                vector = 235,
                matrix = 18,
                array = 88
            ),
            levels = xterm.pal()$paired[c(TRUE, FALSE)],
            range = xterm.pal()$long,
            logical = xterm.pal()$Accent[2:1]
        )
    )
}

#' Automatic styling according to object properties.
#'
#' @param x Object to decide formatting from.
#' @param \dots Sent to \code{\link{style}}.
#' @return See \code{\link{style}}.
#' @seealso style, style.palette
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
style.auto <- function(x, ...){
    pal <- getOption("style.palette", get.palette("dark"))
    
    fg <- pal$fg[
        if(is.null(x)) "null"
        else if(is.character(x)) "character"
        else if(is.numeric(x)) "numeric"
        else if(is.factor(x)) "factor"
        else if(is.logical(x)) "logical"
        else if(is.list(x)) "list"
        else if(is.function(x)) "function"
        else ""
    ]

    ld <- max(1, min(3, length(dim(x))))
    bg <- pal$bg[switch(ld, if(length(x) > 1) "vector" else "scalar", "matrix", "array")]

    style(..., fg=fg, bg=bg)
}


