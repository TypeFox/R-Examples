#' Continuous color palettes.
#' 
#' Functions generating color palettes useable with \R graphical functions, in particular with 
#' \code{\link[sp]{spplot}}. These palettes are continuous, contrast being accentuated by darkening 
#' and lightening extrem values. \code{prevR.demo.pal} plot the available palettes. 
#' \code{prevR.colors.qgis.pal} export a palette in a text file readable by Quantum GIS, 
#' an open-source mapping software.
#' 
#' @param n number of different colors in the palette.
#' 
#' @details \code{\link{prevR.colors.red}} produces a color gradation from white/yellow to red/dark red.\cr
#' \code{\link{prevR.colors.blue}} produces a color gradation from light blue to dark blue.\cr
#' \code{\link{prevR.colors.green}} produces a color gradation from light green to dark green.\cr
#' \code{\link{prevR.colors.gray}} produces a color gradation from white/light gray to dark gray/black.\cr
#' 
#' Functions with a suffix \emph{.inverse} produce the same color gradation, but from dark colors to light ones.
#' 
#' @return 
#' \code{\link{prevR.demo.pal}} plot the color palettes.
#' 
#' \code{\link{prevR.colors.qgis.pal}} export a color palette in a texte file readable by Quantum GIS.
#' 
#' The other functions return a list of colors coded in hexadecimal.
#' 
#' @note 
#' To obtain the liste of colors in RGB (Red/Green/Blue), use the function 
#' \code{\link[grDevices]{col2rgb}}\{\pkg{grDevices}\}. 
#' The code of \code{\link{prevR.demo.pal}} was adapted from the function \code{demo.pal} 
#' presented in the examples of \code{\link[grDevices]{rainbow}}.
#' 
#' @seealso Other color palettes are available in  \R. See for example 
#' \code{\link[grDevices]{rainbow}}\{\pkg{grDevices}\} or the package \pkg{RColorBrewer}.
#' 
#' @examples 
#' prevR.demo.pal(25)
#' prevR.colors.red(5)
#' col2rgb(prevR.colors.red(5))
#' 
#' \dontrun{
#'  prevR.colors.qgis.pal('palette.txt', seq(0,25,length.out=100), 'red')
#' }
#'
#' @name prevR.colors
#' @export
#' @keywords color
prevR.colors.blue <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%2.75
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 26/48, s = seq(to = 1 - 1/(2 * j), 
            from = 1/(2 * j), length = j), v = 1), hsv(h = seq(26/48, 
            33/48, length = i)), if (j > 0) hsv(h = 33/48, v = seq(from = 1 - 
            1/(2 * j), to = 1/(2 * j), length = j), s = 1))
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.blue.inverse <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%2.75
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 33/48, v = seq(to = 1 - 1/(2 * j), 
            from = 1/(2 * j), length = j), s = 1), hsv(h = seq(33/48, 
            26/48, length = i)), if (j > 0) hsv(h = 26/48, s = seq(from = 1 - 
            1/(2 * j), to = 1/(2 * j), length = j), v = 1))
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.gray <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        i <- n - 1
        gray(i:0/i)
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.gray.inverse <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        i <- n - 1
        gray(0:i/i)
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.green <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%2
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 1/3, s = seq(to = 1 - 1/(2 * j), 
            from = 1/(2 * j), length = j), v = 1), hsv(h = seq(1/3, 
            1/3, length = i)), if (j > 0) hsv(h = 1/3, v = seq(from = 1 - 
            1/(2 * j), to = 1/(2 * j), length = j), s = 1))
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.green.inverse <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%2
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 1/3, v = seq(to = 1 - 1/(2 * j), 
            from = 1/(2 * j), length = j), s = 1), hsv(h = seq(1/3, 
            1/3, length = i)), if (j > 0) hsv(h = 1/3, s = seq(from = 1 - 
            1/(2 * j), to = 1/(2 * j), length = j), v = 1))
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.red <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%3.5
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 1/6, s = seq(to = 1 - 1/(2 * j), 
            from = 1/(2 * j), length = j), v = 1), hsv(h = seq(1/6, 
            0, length = i)), if (j > 0) hsv(h = 0, v = seq(from = 1 - 
            1/(2 * j), to = 1/(2 * j), length = j), s = 1))
    }
    else character(0)
}

#' @export
#' @rdname prevR.colors
prevR.colors.red.inverse <-
function (n) 
{
    if ((n <- as.integer(n[1])) > 0) {
        j <- n%/%3.5
        i <- n - 2 * j
        c(if (j > 0) hsv(h = 0, v = seq(to = 1 - 1/(2 * j), from = 1/(2 * 
            j), length = j), s = 1), hsv(h = seq(0, 1/6, length = i)), 
            if (j > 0) hsv(h = 1/6, s = seq(from = 1 - 1/(2 * 
                j), to = 1/(2 * j), length = j), v = 1))
    }
    else character(0)
}

# Cette fonction permet d'afficher un graphique presentant les differentes palettes de prevR
# Exemple : prevR.demo.pal(15)

#' @param border border color.
#' @param main title.
#' @export
#' @rdname prevR.colors
prevR.demo.pal <-
function(n, border = if (n<32) "light gray" else NA, main = NULL) {
    #Fonction basee sur demo.pal - cf. aide de la fonction rainbow : ?rainbow
    ch.col = c(
      "prevR.colors.red.inverse(n)","prevR.colors.blue.inverse(n)", 
      "prevR.colors.green.inverse(n)","prevR.colors.gray.inverse(n)","prevR.colors.red(n)",
      "prevR.colors.blue(n)", "prevR.colors.green(n)", "prevR.colors.gray(n)")
    nt <- length(ch.col)
    i <- 1:n; j <- n / nt; d <- j/6; dy <- 2*d
    if(missing(main)) main= gettextf("prevR.colors palettes for n=%d",n,domain="R-prevR")
    plot(i,i+d, type="n", yaxt="n", xlab="n", ylab="", main=main)
    for (k in 1:nt) {
        rect(i-.5, (k-1)*j+ dy, i+.4, k*j,col = eval(parse(text=ch.col[k])), border = border)
        text(2*j,  k * j +dy/4, ch.col[k],  cex=0.8)
    }
}
# Fonction prevR.colors.qgis.pal() 
# Permet de generer une palette de couleurs utilisable par Quantum GIS
# Il s'agit d'un fichier txt
# Les lignes sont de la forme :
# -12.5,0,68,27,255,Color entry 1
# valeur,Rouge,Vert,Bleu,Transparence,Nom couleur
# pal vaut red, green, blue ou gray (selection de la palette)
# inverse : si true, prendre l'ordre inverse
# at : liste des valeurs de la palette
# n : nombre de couleurs
# file : nom du fichier
# Exemple d'utilisation : prevR.colors.qgis.pal('palette.txt',seq(0,25,length.out=100),'red')

#' @param file file name with extension.
#' @param at list of values of the palette.
#' @param pal color palette to use ("red", "green", "blue" or "gray").
#' @param inverse use the inverse palette?
#' @export
#' @rdname prevR.colors
prevR.colors.qgis.pal <- function (file, at, pal="red", inverse=FALSE) {
    n.pal <- paste('prevR.colors.',pal,sep='')
    if (inverse)
        n.pal <- paste(n.pal,'.inverse',sep='')
    f.pal <- get(n.pal)
    write.table(
        format(data.frame(
            at,
            t(col2rgb(f.pal(length(at)))),
            rep(255,length.out=length(at)),
            rep(n.pal,length.out=length(at))
        ),scientific=FALSE, trim=TRUE),
        file,
        sep=',',dec='.',
        col.names=FALSE, row.names=FALSE, quote=FALSE
    )
}