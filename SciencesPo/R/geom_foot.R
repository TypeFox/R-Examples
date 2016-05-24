#' @encoding UTF-8
#' @title Add Footnote to a ggplot Object
#'
#' @description Add footnotes to \pkg{ggplot2} objects.
#'
#' @param text any text or empty to use default.
#' @param fontsize the font size \code{text}.
#' @param rotn the rotation for the footnote, default is \code{rotation=90}.
#' @param color the color for \code{text}.
#' @param just the justification method.
#'
#' @details At this stage, this function only works for a ggplot object.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @examples
#' # setup data
#' set.seed(51)
#' supply <- rnorm(100,mean=15-seq(1,6,by=.05),sd=1)
#' demand <- rnorm(100,mean=4+seq(1,21,by=.2),sd=.5)
#' time<-seq(1,100,by=1)
#' data <- data.frame(time, supply,demand)
#'
#' # make the plot
#' library(ggplot2)
#' ggplot(data,aes(time)) +
#' geom_line(aes(y=demand),size=1.6, color="#008fd5") +
#' geom_line(aes(y=supply),size=1.6, color="#ff2700") +
#' theme_fte() +
#' annotate("text",x=90,y=12,label="Demand") +
#' annotate("text",x=80,y=23,label="Supply")
#' geom_foot("danielmarcelino.github.io", color = "#77ab43", rotn = -90, just ="right" )
#'
#' @keywords Graphs
#'
#' @import ggplot2
#'
#' @export
#'
`geom_foot` <-
  function(text=NULL, fontsize=10, color=NULL, rotn = 0, just = c("right", "bottom")) {
    if(!is.null(text)){
      text = paste(text)
    } else{
      text = paste(Sys.info()["user"],
                   format(Sys.time(), "%d %b %Y"), sep = " " )
    }
    if(is.null(fontsize)){
      fontsize = .75
    }
    if(is.null(color)){
      color = grDevices::grey(.65)
    }
    #grid::pushViewport(grid::viewport())
    grid::grid.text(label = text ,
              x = 0.99,
              y = 0.01,
              just = just,
              rot = rotn,
              gp = grid::gpar(fontsize = fontsize, col = color))
   # grid::popViewport()
  }
NULL



#' @encoding UTF-8
#' @title Add Title and Subtitle to a ggplot Object
#'
#' @description A production function to make it easy to add title and subtitle to \pkg{ggplot2} objects.
#'
#' @param title A character string as title.
#' @param subtitle A character string as subtitle.
#' @export
#' @examples
#' # setup data
#' set.seed(51)
#' supply <- rnorm(100,mean=15-seq(1,6,by=.05),sd=1)
#' demand <- rnorm(100,mean=4+seq(1,21,by=.2),sd=.5)
#' time<-seq(1,100,by=1)
#' data <- data.frame(time, supply,demand)
#'
#' # make the plot
#' library(ggplot2)
#' ggplot(data,aes(time)) +
#' geom_line(aes(y=demand),size=1.6) +
#' geom_line(aes(y=supply),size=1.6) +
#' annotate("text",x=90,y=12,label="Demand",colour="red") +
#' annotate("text",x=80,y=23,label="Supply",colour="blue") +
#' plotTitleSubtitle("My Title", "My Subtitle")
#'
`plotTitleSubtitle` <- function(title, subtitle = "") {
  ggplot2::ggtitle(bquote(atop(.(title), atop(.(subtitle)))))
}
NULL
