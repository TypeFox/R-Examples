#' ICES Areas
#' 
#' Draw a map showing ICES areas 1--14.
#' 
#' 
#' @param labels is how to annotate areas: \code{"roman"} (I--XIV),
#' \code{FALSE} (no labels), vector (custom), otherwise 1--14.
#' @param diagrams is whether to show diagrams on console.
#' @param col is passed to \code{map}.
#' @param lwd is passed to \code{lines}.
#' @param col.lines is passed to \code{lines}.
#' @param font is passed to \code{text}.
#' @param col.text is passed to \code{text}.
#' @param cex is passed to \code{text}.
#' @return Invisible data frame with coordinates.
#' @note The coordinates were inferred from official ICES maps.
#' @author Arni Magnusson.
#' @keywords hplot spatial utilities
#' @examples
#' 
#' ices()
#' 
#' @export ices
ices <-
function(labels="roman", diagrams=FALSE, col="black", lwd=3, col.lines="grey80", font=2, col.text="orangered", cex=1)
{
  areaLines <- function(area, ...)
  {
    Lines <- function(x) lines(x$East, x$North, ...)
    lapply(area, Lines)
    invisible(NULL)
  }

  coords <- data.frame(Area=as.integer(rep.int(c(1:10,12,14),c(6,9,5,10,11,9,8,4,4,5,9,8))),
                       Line=rep.int(c("0","2","1","4-5-14","4","2-5-6","3","7","all","5-7-12","7","4","6","6-8-12",
                         "all","0","2-5-12"),c(2,4,4,5,5,3,5,2,11,7,2,2,2,4,22,2,6)),
                       North=c(90,68.5,90,72,72,71.2,90,72,72,71.2,90,63,63,62,62,58,57.5,57.5,57,57,58.6,62,62,58,57.5,
                         57.5,57,57,51,51,68,68,63,63,60.5,60.5,60,60,62,62,68,54.5,54.5,60,60,60.5,60.5,58.6,55,55,51,
                         51,55,55,54.5,54.5,48,48,48,48,43,43,43,43,36,36,48,48,36,36,48,62,62,60,60,48,48,59,59,62,90,
                         83.4,90,68,68,59,59,59.8),
                       East=c(51,51,30,30,26,26,30,30,26,26,-11,-11,-4,-4,5,7,7,8,8,8.4,-4,-4,5,7,7,8,8,8.4,1,2,-27,-11,
                         -11,-4,-4,-5,-5,-15,-15,-27,-27,-8.3,-18,-18,-5,-5,-4,-4,-5.9,-5.2,1,2,-5.9,-5.2,-8.3,-18,-18,
                         -4.6,-4.6,-18,-18,-9.3,-9.3,-18,-18,-5.6,-42,-18,-18,-42,-42,-27,-15,-15,-18,-18,-42,-42,-27,
                         -27,-40,-40,-11,-11,-27,-27,-44,-44))
  chunks <- lapply(split(coords,coords$Area), function(x) split(x,x$Line))
  map(xlim=c(-55,55), ylim=c(30,90), fill=TRUE, col=col)
  sapply(1:12, function(i) areaLines(chunks[[i]], lwd=lwd, col=col.lines))
  if(!identical(labels, FALSE))
  { #             I    II   III    IV      V     VI    VII   VIII     IX    X  XII   XIV
    textN <- c(72.4, 72.4, 56.4, 56.4,  63.2,  57.3,  51.3,  45.5,  39.7,  43 , 54, 72.4)
    textE <- c(  40,    5,   19,    3, -13.6, -13.6, -13.6, -13.6, -13.6, -32, -32,  -32)
    strings <-
      if(identical(labels,"roman"))
        as.roman(c(1:10,12,14))
      else if(length(labels) > 1)
        labels
      else
        c(1:10,12,14)
    text(textE, textN, strings, font=font, col=col.text, cex=cex)
  }
  if(diagrams)
  {
    cat("",
        "Area I: Barents Sea",
        "",
        "             + 90N,30E   + 90N,51E",
        "             |           |",
        "             |           |",
        "  72N,26E +--+ 72N,30E   |",
        "          |              |",
        "          |              |",
        "71.2N,26E +              |",
        "                         |",
        "                         |",
        "                         + 68.5N,51E",
        "",
        "",
        "",
        "Area II: Norwegian Sea",
        "",
        "90N,11W +                                         + 90N,30E",
        "        |                                         |",
        "        |                                         |",
        "        |                              72N,26E +--+ 72N,30E",
        "        |                                      |",
        "        |                                      |",
        "        |                            71.2N,26E +",
        "        |",
        "        |",
        "63N,11W +--------+ 63N,4W",
        "                 |",
        "                 |",
        "          62N,4W +--------+ 62N,5E",
        "",
        "",
        "",
        "Area III: Baltic",
        "",
        "  58N,7E +",
        "         |",
        "         |",
        "57.5N,7E +-----------+ 57.5N,8E",
        "                     |",
        "                     |",
        "              57N,8E +----------+ 57N,8.4E",
        "",
        "",
        "",
        "Area IV: North Sea",
        "",
        "  62N,4W +--------------------+ 62N,5E",
        "         |",
        "         |",
        "58.6N,4W +",
        "                                           58N,7E +",
        "                                                  |",
        "                                                  |",
        "                                         57.5N,7E +--------+ 57.5N,8E",
        "                                                           |",
        "                                                           |",
        "                                                    57N,8E +----------+ 57N,8.4E",
        "           51N,1E +--+ 51N,2E",
        "",
        "",
        "",
        "Area V: Iceland and Faroes",
        "",
        "68N,27W +-----------------------------+ 68N,11W",
        "        |                             |",
        "        |                             |",
        "        |         .    .    . 63N,11W +-----------------------------+ 63N,4W",
        "        |                                                           |",
        "        |         .                                                 |",
        "62N,27W +---------+ 62N,15W                                         |",
        "                  |                                                 |",
        "                  |                                                 |",
        "                  |                               60.5N,5W +--------+ 60.5N,4W",
        "                  |                                        |",
        "                  |                                        |",
        "          60N,15W +----------------------------------------+ 60N,5W",
        "",
        "",
        "",
        "Area VI: Scotland",
        "",
        "                                                     60.5N,5W +----------+ 60.5N,4W",
        "                                                              |          |",
        "                                                              |          |",
        "  60N,18W +---------------------------------------------------+ 60N,5W   |",
        "          |                                                              |",
        "          |                                                              |",
        "          |                                                              + 58.6N,4W",
        "          |",
        "          |",
        "          |                 55N,5.9W ---- 55N,5.2W",
        "          |",
        "          |",
        "54.5N,18W +--+ 54.5N,8.3W",
        "",
        "",
        "",
        "Area VII: Sole",
        "",
        "                            55N,5.9W +--+ 55N,5.2W",
        "",
        "54.5N,18W +--+ 54.5N,8.3W",
        "          |",
        "          |",
        "          |",
        "          |",
        "          |                                                     51N,1E +--+ 51N,2E",
        "          |",
        "          |",
        "          |",
        "  48N,18W +----------------------------------------+ 48N,4.6W",
        "",
        "",
        "",
        "Area VIII: Biscay",
        "",
        "48N,18W +-------------+ 48N,4.6W",
        "        |",
        "        |",
        "43N,18W +--+ 43N,9.3W",
        "",
        "",
        "",
        "Area IX: Portugal",
        "",
        "43N,18W +--+ 43N,9.3W",
        "        |",
        "        |",
        "36N,18W +-------------+ 36N,5.6W",
        "",
        "",
        "",
        "Area X: Azores",
        "",
        "48N,42W +--+  48N,18W",
        "        |  |",
        "        |  |",
        "36N,42W +--+  36N,18W",
        "",
        "",
        "",
        "Area XII: North Azores",
        "",
        "          62N,27W +-----------------------------+ 62N,15W",
        "                  |                             |",
        "                  |                             |",
        "                  |           60N,18W +---------+ 60N,15W",
        "                  |                   |",
        "                  |                   |",
        "59N,42W +---------+ 59N,27W           |",
        "        |                             |",
        "        |                             |",
        "48N,42W +-----------------------------+ 48N,18W",
        "",
        "",
        "",
        "Area XIV: Greenland",
        "",
        "              90N,40W +                   + 90N,11W",
        "                      |                   |",
        "                      |                   |",
        "            83.4N,40W +                   |",
        "                                          |",
        "                                          |",
        "                        68N,27W +---------+ 68N,11W",
        "                                |",
        "                                |",
        "59.8N,44W +                     |",
        "          |                     |",
        "          |                     |",
        "  59N,44W +---------------------+ 59N,27W",
        "",
        "",
        "",
        " 1 Barents Sea",
        " 2 Norwegian Sea",
        " 3 Baltic",
        " 4 North Sea",
        " 5 Iceland and Faroes",
        " 6 Scotland",
        " 7 Sole",
        " 8 Biscay",
        " 9 Portugal",
        "10 Azores",
        "12 North Azores",
        "14 Greenland",
        sep="\n")
  }
  invisible(coords)
}
