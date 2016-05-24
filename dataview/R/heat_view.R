#' Display heatmaps and heatvectors.
#'
#' Quickly see the overall pattern of a variable in the terminal.
#'
#' @param x Vector to be displayed.
#' @param pal Palette. Either the name of a palette defined in \code{\link[xtermStyle]{xterm.pal}}
#'   or an integer vector with color indices (see \code{\link[xtermStyle]{display.xterm.colors}}).
#' @param range The numerical range which the palette describes. See \code{\link[xtermStyle]{discrete.color}}
#'   for more info.
#' @param mark Single letter marks to be displayed on top of the color.
#' @return Nothing
#' @examples
#' data(iris)
#' heat.view(iris$Species)
#' heat.view(matrix(iris$Petal.Width, 3, 50, byrow=TRUE,
#'                  dimnames=list(levels(iris$Species), NULL)), pal="purples")
#'
#' run.status <- factor(runif(100) < .95, labels=c("Fail", "Pass"))
#' heat.view(run.status, pal=1:2)
#'
#' #Tip for displayig the element names of a named vector:
#' a <- runif(7)
#' names(a) <- c("ATM", "CHK1", "CDC25", "p53", "CDC2", "CDK2", "CDK4")
#' heat.view(a)            # No names displayed
#' heat.view(as.matrix(a)) # Names displayed
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
heat.view <- function(x, pal, range, mark=NULL){
    terminal.width <- if(Sys.getenv("COLUMNS") == ""){
        getOption("width", 80L)
    } else {
        as.integer(Sys.getenv("COLUMNS"))
    }
    if(inherits(x, "Surv")){
        if(missing(mark)) mark <- ifelse(x[,2] == 0, " ", as.character(x[,2]))
        x <- x[,1]
    }
    n <- length(x)
    nc <- nchar(n)
    row.length <- floor((terminal.width - nc - 2)/10)*10
    if(is.null(mark))
        mark <- ifelse(is.na(x), "-", " ")
    mark <- substr(paste0(mark, " "), 1, 1)

    switch(class(x),
        `logical` = {
            if(missing(pal)) pal <- getOption("style.palette")$logical
            legend <- mapply(style, c("FALSE", "TRUE"), fg=pal)
            legend <- paste(legend, collapse=", ")
            x <- x + 1
        },
        `factor` = {
            if(missing(pal)) pal <- getOption("style.palette")$levels
            if(length(levels(x)) > length(pal)){
                legend <- mapply(style, 
                    c(head(levels(x), length(pal)-2), "...", tail(levels(x), 1)),
                    fg = pal)
                x <- findInterval(as.integer(x), c(1:(length(pal)-1), length(levels(x))))
            } else {
                legend <- mapply(style, levels(x), fg=head(pal, length(levels(x))))
                x <- as.integer(x)
            }
            legend <- paste(legend, collapse=", ")
        }, {
            if(missing(pal)) pal <- getOption("style.palette")$range
            if(missing(range)) range <- c(-1, 1)*max(abs(x))
            cuts <- seq(range[1], range[2], length.out=length(pal)+1)
            x <- findInterval(x, cuts[-(length(pal)+1)])
            legend <- c(sprintf("%g ", cuts[min(x, na.rm=TRUE)]),
                        mapply(style, "=", fg=pal[seq(min(x), max(x))]),
                        sprintf(" %g", cuts[max(x, na.rm=TRUE)]))
            legend <- paste(legend, collapse="")
        }
    )
    mark <- split(mark, gl(ceiling(n/row.length), row.length, n))
    x <- split(x, gl(ceiling(n/row.length), row.length, n))
    cat(sep="\n",
        sprintf("%s: %s%s",
            sprintf(sprintf("%%%ii", nc), 1 + (seq_along(x)-1)*row.length),
            mapply(function(x, m) paste(mapply(style, m, bg=pal[x]), collapse=""), x, mark),
            style.clear()
        )
    )
    cat("\n", legend, "\n")
}


