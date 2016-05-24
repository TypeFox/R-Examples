#' Display color tables
#'
#' The xterm colour table consist of the ANSI colours (16), the web
#' colour cube (216) and additional shades of grey not including full white and
#' black (16). However these are not strictly defined but can vary somewhat
#' between systems and configurations.
#'
#' @param numbers Logical, whether to display colour indices.
#' @param perm Rotation of the colour cube, supplied as a permutation of its
#'   dimensions. Sent to \code{\link[base]{aperm}}.
#' @return Nothing
#' @examples
#' display.xterm.colors()
#' display.ANSI.colors()
#' @seealso style
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
display.xterm.colors <- function(perm=1:3) {
    my.mode <- getOption("style.mode")
    tryCatch({
        ## 16 ANSI colours
        for(i in 0:1) {
            for(j in 0:7){
                cat(style(fg = 8 * as.numeric(i == 0 & j == 0), bg = i * 8 + j, mode="xterm-256color"))
                cat(sprintf("%4i", i*8+j))
            }
            style.clear()
            cat("\n")
        }
        cat("\n")
        
        ## 216 Web colours
        color.cube <- aperm(array(16:231, c(6,6,6)), perm=perm)
        for(i in 1:6^3){
            x <- (i-1) %% 6 + 1
            y <- floor((i-1)%%72 / 12) + 1
            z <- 2*floor((i-1)/72) + ((i-1)%%12 >= 6) + 1
            my.color = color.cube[x,y,z]

            cat(style(fg = 8 * as.numeric(my.color == 16), bg = my.color, mode="xterm-256color"))
            cat(sprintf("%4s", my.color))
            if(i %% 6 == 0 && i %% 12 != 0)
                cat(style.clear(), "  ")
            if(i %% (2*6) == 0)
                cat(style.clear(), "\n")
            if(i %% (2*6*6) == 0)
                cat("\n")
        }
        
        ## 16 
        for(i in 0:1) {
            for(j in 0:11){
                cat(style(fg = (238 + j) * as.numeric(i == 0), bg = 232 + i * 12 +j, mode="xterm-256color"))
                cat(sprintf("%4s", 232 + i * 12 +j))
            }
            cat(style.clear(), "\n")
        }
        options(style.mode = my.mode)
    }, interrupt = function(...) cat(style.clear()))
}
#' @rdname display.xterm.colors
#' @export
display.ANSI.colors <- function(numbers=TRUE){
    tryCatch({
        my.mode <- getOption("style.mode")
        ## 16 ANSI colours
        for(i in 0:1) {
            for(j in 0:7){
                cat(style(fg = 8 * as.numeric(i == 0 & j == 0), bg = i * 8 + j, mode="ANSI"))
                cat(sprintf("%4s", i * 8 + j))
            }
            style.clear()
            cat("\n")
        }
        cat("\n")
        options(style.mode = my.mode)
    }, interrupt = function(...) cat(style.clear()))
}


#' Get predefined colour palettes
#'
#' All except "GnRd" and "long" are basen on the color brewer palettes, see
#' \code{\link[RColorBrewer]{brewer.pal}} of the \code{RColorBrewer} package.
#'
#' @return A list of vectors with colour indices.
#' @examples
#' display.xterm.pal()
#' pal <- xterm.pal()$Accent
#'
#' freqs <- runif(6)
#' fruits <- factor(sample(6, size=30, replace=TRUE, freqs/sum(freqs)),
#'                  labels=c("apple", "grapes", "banana", "lemon",
#'                           "blueberry", "raspberry"))
#' for(i in 1:length(fruits))
#'     cat(style(fruits[i], "\n", fg=pal[fruits[i]]))
#'
#' @seealso display.xterm.pal, display.xterm.colors
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
xterm.pal <- function(){
    list(
        YlOrRd = c(230, 228, 221, 208, 202, 196, 160, 124, 88),
        YlOrBr = c(230, 229, 228, 222, 215, 209, 166, 130, 94),
        YlGnBu = c(230, 194, 156, 115, 80, 38, 26, 20, 18),
        YlGnBu = c(230, 194, 156, 115, 80, 38, 26, 20, 18),
        YlGn = c(230, 228, 156, 113, 77, 41, 35, 28, 22),
        Reds = c(231, 224, 217, 211, 203, 196, 124, 88, 52),
        RdPu = c(225, 219, 213, 207, 201, 164, 127, 90, 53),
        Purples = c(231, 189, 183, 147, 141, 105, 99, 63, 57),
        PuRd = c(225, 183, 177, 170, 201, 163, 126, 89, 52),
        PuBuGn = c(159, 117, 74, 39, 38, 37, 36, 29, 23),
        PuBu = c(159, 117, 81, 75, 38, 32, 26, 21, 18),
        OrRd = c(230, 222, 215, 209, 203, 196, 160, 124, 88),
        Oranges = c(230, 222, 215, 214, 208, 202, 166, 130, 94),
        Greys = c(254, 251, 249, 246, 243, 240, 237, 235, 234),
        Greens = c(194, 157, 120, 84, 77, 71, 34, 28, 22),
        GnBu = c(194, 157, 84, 42, 37, 31, 26, 20, 18),
        BuPu = c(159, 117, 75, 69, 63, 57, 56, 55, 53),
        BuGn = c(159, 117, 75, 39, 37, 36, 35, 28, 22),
        Blues = c(159, 117, 75, 39, 33, 27, 21, 19, 17),
        
        Set3 = c(80, 229, 61, 203, 75, 208, 76, 218, 248, 93, 158, 220),
        Set2 = c(79, 209, 104, 182, 155, 221, 180, 245),
        Set1 = c(160, 27, 35, 93, 208, 227, 130, 213, 244),
        pastel2 = c(122, 223, 117, 225, 158, 228, 224, 250),
        pastel1 = c(217, 117, 193, 147, 223, 230, 187, 225, 254),
        paired = c(117, 33, 157, 78, 218, 160, 216, 202, 183, 93, 229, 130),
        Dark2 = c(30, 130, 62, 162, 70, 179, 94, 243),
        Accent = c(78, 141, 221, 229, 26, 162, 130, 243),
            
        Spectral = c(53, 125, 203, 216, 222, 228, 191, 114, 36, 25, 57),
        RdYlGn = c(124, 160, 202, 214, 221, 228, 191, 118, 76, 34, 22),
        RdYlBu = c(124, 160, 202, 214, 221, 228, 159, 117, 75, 33, 26),
        RdGy = c(88, 124, 196, 203, 217, 231, 251, 246, 242, 237, 234),
        RdBu = c(88, 124, 196, 203, 217, 255, 159, 117, 75, 33, 26),
        PuOr = c(94, 166, 209, 215, 222, 255, 225, 177, 135, 92, 54),
        PRGn = c(54, 92, 135, 177, 189, 255, 194, 119, 76, 34, 22),
        PiYG = c(164, 171, 213, 219, 225, 255, 194, 119, 76, 34, 22),
        BrBG = c(94, 166, 209, 215, 222, 255, 158, 115, 73, 30, 23),
            
        GnRd = c(seq(46,226,36), seq(220,196,-6)),
        long = c(18:20, seq(21,51,by=6), seq(87,195,by=36), 231:226, seq(220,196,by=-6), seq(160,88,by=-36)),
        DownUp = c(seq(51, 21, by=-6), 20:16, seq(52, 196, by=36), seq(202, 226, by=6)),
        BuPuYl = c(seq(87, 57, by=-6), 56:53, seq(89, 197, by=36), seq(203, 227, by=6)),
        jet = c(17:21, 27,33,39,45, 51:46, 82,118,154,190,226, 220,214,208,202, 196:201, 207,213,219,225,231)
    )
}
#' @rdname xterm.pal
#' @export
xterm.pal.inv <- function(){
    lapply(xterm.pal(), function(p){
        new.r <- 5 - floor((p-16)/36)
        new.g <- 5 - floor((p-16) %% 36 / 6)
        new.b <- 5 - (p-16) %% 6
        new.ANSI <- 15 - p
        new.grey <- 255 - p + 232
        p[p > 15 & p < 232] <- (new.r * 36 + new.g * 6 + new.b + 16)[p > 15 & p < 232]
        p[p < 16] <- new.ANSI[p < 16]
        p[p > 231] <- new.grey[p > 231]
        p
    })
}
#' @rdname xterm.pal
#' @export
display.xterm.pal <- function(){
    pal <- xterm.pal()

    nc <- max(nchar(names(pal)))
    terminal.width <- if(Sys.getenv("COLUMNS") == ""){
        getOption("width", 80)
    } else {
        as.integer(Sys.getenv("COLUMNS"))
    }
    cols.per.line <- floor((terminal.width - nc - 2) / 4)
    tryCatch({
        for(i in 1:length(pal)){
            cat(sprintf(sprintf("\n%%%is: ", nc), names(pal)[i]))
            for(j in 1:length(pal[[i]])){
                if(j > 1 && j %% cols.per.line == 1)
                    cat(style.clear(), sprintf(sprintf("\n%%%is", nc+2), ""))
                cat(style(bg=pal[[i]][j]), sprintf("%4i", pal[[i]][j]), sep="")
            }
            cat(style.clear(), "\n", sep="")
        }
    }, interrupt = function(...) cat(style.clear()))
}

#' Map numbers onto a palette
#' 
#' The continuous interval defined by \code{range} is divided into bins of
#' equal size. Each bin is mapped to a colour in the palette defined by \code{pal}.
#' The values in \code{x} are then assigned to the bins and their corresponding
#' colours are returned. Values outside the interval are assigned to the border
#' bins.
#' 
#' @param x Continuous numbers.
#' @param range The interval in \code{x} that will be mapped to the palette.
#' @param pal Palette. Can be the name of a predefined palette, as returned by
#'   \code{\link{xterm.pal}}, or a vector of colour indices directly.
#' @return Colour indices from \code{pal} corresponding to where in the range
#'   the values in \code{x} are.
#' @examples
#' error.rates <- .6*runif(10)
#' for(q in error.rates)
#'   style(q, "\n", fg=discrete.color(q, c(0, .5), "GnRd"))
#' @seealso xterm.pal
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
discrete.color <- function(x, range=range(x), pal="GnRd"){
    if(is.character(pal)) pal <- xterm.pal()[[pal]]
    n.pal <- length(pal)
    pal[findInterval(x, seq(range[1], range[2], length=n.pal+1)[2:n.pal])+1]
}

#' Convert xterm256 color code to ANSI code
#' 
#' @param x Integer specifying xterm256 color.
#' @return Integer that approximates \code{x} in the ANSI palette.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
xterm256.to.ANSI <- function(x){
    x <- as.integer(x)
    if(x < 16){
        x %% 8
    } else if(x < 232){
        binary.color <- round( c((x-16) %% 6, floor((x-16) %% 36 / 6), floor((x-16) / 36)) / 5)
        # Correspons to blue, green, red
        sum(c(4,2,1)*binary.color)
    } else {
        7 * (x > 243)
    }
}
