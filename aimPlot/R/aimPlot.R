##' Pie Like Plot for Completeness
##'
##' Create a pie like plot to visualize how close a project is to completeness or
##' achievement of the aims. The middle point show the total completeness of the
##' aim.  Imagine it's like a dartboard and the center means 100\% completeness/achievement.
##' The standard distribution for the proportion of completeness allocated in the pie is 50\%, 80\% and 100\%.
##' Items for aim should be on the first row and the percentage of each items is on the second row in the data set.
##'
##' @note The \code{ggplot2} package is required to run this function
##' @param data Data set
##' @param title Title for the plot
##' @param size Size of the point
##' @param pct1 Percentage first pie proportion
##' @param pct2 Percentage second pie proportion
##' @param col1 Colour of the first pie proportion
##' @param col2 Colour of the second pie proportion
##' @param col3 Colour of the third pie proportion
##' @details These parameters should be in the dataframe:
##' \itemize{
##'  \item \code{Variable} {1st column: The various aims to be completed}
##'  \item \code{Percentage} {2nd column: The percentage of completeness}
##' }
##' @source  \code{RegData} is example data set from Norwegian Diabetes Registry
##' @author Yusman Kamaleri, \email{ybkamaleri@gmail.com}
##'
##' @examples
##'
##' # basic usage
##' library("aimPlot")
##' aimPlot(data = RegData, title = "Plot title")
##' aimPlot(RegData, "Plot title", size=10, col1="blue", col2="green", col3="yellow")
##' aimPlot(RegData, pct1 = 20, pct2 = 60)
##' 
##' @export

aimPlot <- function(data, title, size, pct1, pct2, col1, col2, col3) {

    if (missing(data)) {
        stop("'data' must be provided",
             call. = FALSE)
    }
    
    if (missing(title)) {
        title <- " "
    }

    if (missing(size)) {
        size <- 14
    }

    if (missing(pct1)) {
        pct1 <- 50
    }
    
    if (missing(pct2)) {
        pct2 <- 80
    }

    if (missing(col2)) {
        col2 <- "#2171B5"
    }

    if (missing(col1)) {
        col1 <- "#99CCFF"
    }

    if (missing(col3)) {
        col3 <- "#000033"
    }
    
    aim <- data[,1]
    percent <- data[,2]
    aim <- gsub(" ", "\n", levels(data[ ,1]))
    pct2a <- 2 + pct2
    
    p <- ggplot2::ggplot(data) +
        ggplot2::scale_x_discrete() +
        ggplot2::scale_y_reverse() +
        ggplot2::geom_rect(xmin=Inf, xmax = -Inf, ymin = 0 - pct2a, ymax = 0 - pct1, fill=col2) +
        ggplot2::geom_rect(xmin=Inf, xmax = -Inf, ymin = 0 - pct1, ymax = 0, fill=col1) +
        ggplot2::geom_rect(xmin=Inf, xmax = -Inf, ymin = -100, ymax = 0 - pct2, fill=col3) +
        ggplot2::geom_vline(xintercept=1:8, size=1.5, color="white") +
        ggplot2::geom_hline(yintercept=c(pct1, pct2, 100), size=0.1, color="white") +
        ggplot2::geom_point(ggplot2::aes(x=aim, y=percent, ymin=0, ymax=100),
                            shape=21, fill="#FF9933", size=size, position="identity") +
        ggplot2::ggtitle(title) +
        ggplot2::coord_polar() +
        ggplot2::theme(
            plot.title=ggplot2::element_text(face = "bold", color = "black", size = 17, ),
            panel.background=ggplot2::element_rect(fill = c("white")),
            panel.grid=ggplot2::element_blank(),
            panel.grid.major=ggplot2::element_line(size=2),
            panel.grid.minor.y=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_text(vjust=5),
            axis.text=ggplot2::element_text(size = 12, color = "black", face = "bold"),
            axis.text.y=ggplot2::element_blank(),
            axis.title.x=ggplot2::element_blank(),
            axis.title.y=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank())

    return(p)
}


