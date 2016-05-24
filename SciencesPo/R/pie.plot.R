if(getRversion() >= "2.15.1") globalVariables(c(".data", "pie"))

#' @encoding UTF-8
#' @title  Creates a pie chart using ggplot2.
#'
#' @description Use pie charts with care. See http://www.edwardtufte.com/bboard/q-and-a-fetch-msg?msg_id=00018S
#' on Edward Tufte's website for good arguments against the use of pie charts.
#' For a contrary point-of-view, see Spence's article, No Humble Pie: The Origins and
#' Usage of a Statistical Chart (http://www.psych.utoronto.ca/users/spence/Spence%202005.pdf).
#'
#' @param .data the data frame.
#' @param var the name of the column to generate the pie chart for.
#' @param label The label for the legend.
#' @import ggplot2
#' @examples
#' if (interactive()) {
#' x = sample(10, 100, rep = TRUE)
#' z = sample(letters[1:3],100, rep=TRUE)
#' dat = data.frame(x,z)
#' pie.plot(dat, 'x', 'z')
#' }
#' @export
`pie.plot` <- function(.data, var, label=var) {
 # res <- dplyr::count(.data, ...)
  .data$pie = .data[,var]
  l = levels(.data$pie)
  t = table(.data$pie)
  levels(.data$pie) = paste(names(t), " ", format(100*t/sum(t),digits=1), "%", sep="")
  p = ggplot(.data, aes(x=factor(1), fill=pie)) +
    geom_bar(width=1) +
    coord_polar(theta="y") +
    xlab("") + ylab("") +
    scale_fill_hue(name=label, breaks=levels(.data$pie), labels=levels(.data$pie)) + theme_pub() +
    theme(axis.text=element_blank(), axis.ticks=element_blank())
  print(p)
}
