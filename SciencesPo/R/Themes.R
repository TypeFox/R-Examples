#' @encoding UTF-8
#' @title The Default Theme
#'
#' @description After loading the SciencesPo package, this theme will be
#' set to default for all subsequent graphs made with ggplot2.
#'
#' @param legend Enables to set legend position, default is "bottom".
#' @param font_family Default font family.
#' @param font_size Overall font size. Default is 14.
#' @param line_width Default line size.
#' @param axis.line.x Enables to set x axis line.
#' @param axis.line.y Enables to set y axis line.
#' @return The theme.
#' @seealso \code{\link[ggplot2]{theme}}, \code{\link{theme_538}}, \code{\link{theme_blank}}.
#' @examples
#' ggplot(diamonds,aes(cut, group=1)) + geom_bar()+
#' geom_freqpoly(stat="count",size=2) + scale_color_pub() + theme_pub(line_width=1)
#'
#' dat <- data.frame()
#' for(i in 1:4)
#' dat <- rbind(dat, data.frame(set=i, x=anscombe[,i], y=anscombe[,i+4]))
#'
#' ggplot(dat, aes(x, y)) + geom_point(size=5, color="red",
#' fill="orange", shape=21) + geom_smooth(method="lm", fill=NA,
#' fullrange=TRUE) + facet_wrap(~set, ncol=2)
#'
#' @export
`theme_pub` <- function(legend = 'bottom',
                        font_family = 'sans',
                        font_size = 13,
                        line_width = .5,
                        axis.line.x = element_line(),
                        axis.line.y = element_blank()){
half_line <- font_size / 2
small_rel <- 0.857
small_size <- small_rel * font_size
theme_grey(base_size = font_size, base_family = font_family) %+replace%
theme(rect = element_rect(fill = "transparent",
                              colour = NA,
                              color = NA,
                              size = 0,
                              linetype = 0),
          text = element_text(family = font_family,
                              face = "plain",
                              colour = "black",
                              size = font_size,
                              hjust = 0.5,
                              vjust = 0.5,
                              angle = 0,
                              lineheight = .9,
                              margin = ggplot2::margin(),
                              debug = FALSE),
axis.text = element_text(color="black", size = small_size),
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(margin = ggplot2::margin(t = small_size / 4), vjust = 1),
          axis.text.y = element_text(margin = ggplot2::margin(r = small_size / 4), hjust = 1),
          axis.title.x = element_text(
            margin = ggplot2::margin(t = small_size / 2, b = small_size / 4)
          ),
          axis.title.y = element_text(
            angle = 90,
            margin = ggplot2::margin(r = small_size / 2, l = small_size / 4),
          ),
axis.ticks = element_line(color = "#525252", size = line_width),
axis.line = element_line(color = "#525252", size = line_width),
          legend.position = legend,
          # legend.position = c(-0.03, 1.05),
          # legend.justification = c("left", "top"),
          legend.direction = "horizontal",
          legend.key = element_rect(color = NA),
          legend.margin = grid::unit(0, "cm"),
          legend.key.size = grid::unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          legend.text  = element_text(size = rel(small_rel)),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_line(color="#F0F0F0"),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face="bold", size = rel(small_rel)),
strip.background = element_rect(fill = "grey80", color = "grey50", size = 0),
# margins starting with top, right, bottom and left
      plot.margin = grid::unit(c(0.4, 0.2, 0.1, 0.1),"cm"),
          plot.background = element_blank(),
          plot.title = element_text(face = "bold",
                                    size = font_size,
                                    margin = ggplot2::margin(b = half_line))
    )
}
NULL





#' @title Themes for ggplot2 Graphs
#'
#' @description  Theme for plotting  with ggplot2.
#'
#' @param legend Enables to set legend position, default is "none".
#' @param font_family Default font family.
#' @param font_size Overall font size. Default is 13.
#' @param colors Default colors used in the plot in the following order: background, lines, text, and title.
#' @return The theme.
#'
#' @examples
#' qplot(1:10, (1:10)^3) + theme_fte()
#'
#' @export
#' @aliases theme_538
`theme_fte` <- function(legend = 'none', font_size = 12, font_family = 'sans', colors=c('#F0F0F0', '#D0D0D0', '#535353', '#3C3C3C') ){
theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      # Base elements which are not used directly but inherited by others

      line = element_line(color = colors[2], size = 0.75,
                          linetype = 1, lineend = "butt"),
      rect = element_rect(fill = colors[1], color = colors[1],
                          size = 0.5, linetype = 1),
      text = element_text(family = font_family, face = 'bold',
                          color = colors[3], size = font_size,
                          hjust = 0.5, vjust = 0.5, angle = 0,
                          lineheight = 1, margin = ggplot2::margin(1,1,1,1), debug = FALSE),
      # Modified inheritance structure of text element
plot.title = element_text(size = rel(1.5), family = '' ,
                                face = 'bold', hjust = -0.05,
                                vjust = 1.5, color = colors[4],
                          margin = ggplot2::margin(), debug = FALSE),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(),
      # Modified inheritance structure of line element
      axis.ticks =  element_line(),
      panel.grid.major =  element_line(),
      panel.grid.minor =  element_blank(),
      # Modified inheritance structure of rect element
      plot.background =  element_rect(),
      panel.background =  element_rect(),
      strip.background = element_rect(),
      legend.background = element_rect(linetype = 0),
      legend.margin = grid::unit(font_size * 1.1, "points"),
      legend.key = element_rect(linetype = 0),
      legend.key.size = grid::unit(1.1, "lines"),
legend.key.height = NULL,
      legend.key.width = NULL,
legend.text = element_text(size = rel(1.1)),
      legend.text.align = NULL,
legend.title = element_text(size = rel(1), hjust = 0),
      legend.title.align = NULL,
      legend.position = legend,
      legend.direction = NULL,
      legend.box = "vertical",
      legend.justification = "center",
      #legend.key =  element_rect(color = '#D0D0D0'),
      # Modifiying legend.position
      complete = TRUE
    )
}
NULL

#' @export
#' @rdname theme_fte
theme_538 <- theme_fte


#' Create a Completely Empty Theme
#'
#' The theme created by this function shows nothing but the plot panel.
#' @param font_family Default font family.
#' @param font_size Overall font size. Default is 12.
#' @return The theme.
#' @examples
#' # plot with small amount of remaining padding
#' qplot(1:10, (1:10)^2) + theme_blank()
#' # remaining padding removed
#' qplot(1:10, (1:10)^2) + theme_blank() + labs(x = NULL, y = NULL)
#' @export
`theme_blank` <- function(font_size = 12, font_family = ""){
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      rect              = element_rect(fill = "transparent", colour = NA, color = NA, size = 0, linetype = 0),
      line              = element_blank(),
      text              = element_blank(),
      title             = element_blank(),
      # to debug, uncomment next line
      #plot.background   = element_rect(colour = "blue", fill = "cyan"),
      panel.background  = element_blank(),
      axis.ticks.length = grid::unit(0, "cm"),
      legend.position   = "none",
      panel.margin      = grid::unit(c(0, 0, 0, 0), "cm"),
      plot.margin       = grid::unit(c(0, 0, 0, 0), "cm")
    )
}
NULL





#' @title Palette data for the themes used by package
#'
#' @description Data used by the palettes in the package.
#'
#' @format A \code{list}.
#'
themes_data <- {
x <- list()
x$pub <- list()
x$pub$colors <-
  list(
  tableau20=c(
    rgb(31, 119, 180, max = 255),
    rgb(174, 199, 232, max = 255),
    rgb(255, 127, 14, max = 255),
    rgb(255, 187, 120, max = 255),
    rgb(44, 160, 44, max = 255),
    rgb(152, 223, 138, max = 255),
    rgb(214, 39, 40, max = 255),
    rgb(255, 152, 150, max = 255),
    rgb(148, 103, 189, max = 255),
    rgb(197, 176, 213, max = 255),
    rgb(140, 86, 75, max = 255),
    rgb(196, 156, 148, max = 255),
    rgb(227, 119, 194, max = 255),
    rgb(247, 182, 210, max = 255),
    rgb(127, 127, 127, max = 255),
    rgb(199, 199, 199, max = 255),
    rgb(188, 189, 34, max = 255),
    rgb(219, 219, 141, max = 255),
    rgb(23, 190, 207, max = 255),
    rgb(158, 218, 229, max = 255)
    ),
  tableau10medium=c("#729ECE",
      "#FF9E4A",
      "#67BF5C",
      "#ED665D",
      "#AD8BC9",
      "#A8786E",
      "#ED97CA",
      "#A2A2A2",
      "#CDCC5D",
      "#6DCCDA"),
  pub12=c(
    rgb(56, 108, 176, max = 255),
    rgb(253, 180, 98, max = 255),
    rgb(127, 201, 127, max = 255),
    rgb(239, 59, 44, max = 255),
    rgb(102, 37, 6, max = 255),
    rgb(166, 206, 227, max = 255),
    rgb(251, 154, 153, max = 255),
    rgb(152, 78, 163, max = 255),
    rgb(255, 255, 51, max = 255)
    ),
  gray5=c("#60636A",
            "#A5ACAF",
            "#414451",
            "#8F8782",
            "#CFCFCF"),
  trafficlight=c("#B10318",
                   "#DBA13A",
                   "#309343",
                   "#D82526",
                   "#FFC156",
                   "#69B764",
                   "#F26C64",
                   "#FFDD71",
                   "#9FCD99"),
  bluered12=c("#2C69B0",
      "#B5C8E2",
      "#F02720",
      "#FFB6B0",
      "#AC613C",
      "#E9C39B",
      "#6BA3D6",
      "#B5DFFD",
      "#AC8763",
      "#DDC9B4",
      "#BD0A36",
      "#F4737A"),
  purplegray12=c("#7B66D2",
      "#A699E8",
      "#DC5FBD",
      "#FFC0DA",
      "#5F5A41",
      "#B4B19B",
      "#995688",
      "#D898BA",
      "#AB6AD5",
      "#D098EE",
      "#8B7C6E",
      "#DBD4C5"),
  greenorange12=c("#32A251",
      "#ACD98D",
      "#FF7F0F",
      "#FFB977",
      "#3CB7CC",
      "#98D9E4",
      "#B85A0D",
      "#FFD94A",
      "#39737C",
      "#86B4A9",
      "#82853B",
      "#CCC94D"),
  cyclic=c("#1F83B4",
      "#1696AC",
      "#18A188",
      "#29A03C",
      "#54A338",
      "#82A93F",
      "#ADB828",
      "#D8BD35",
      "#FFBD4C",
      "#FFB022",
      "#FF9C0E",
      "#FF810E",
      "#E75727",
      "#D23E4E",
      "#C94D8C",
      "#C04AA7",
      "#B446B3",
      "#9658B1",
      "#8061B4",
      "#6F63BB"),
  colorblind=c(
    rgb(0, 107, 164, max = 255),
    rgb(255, 128, 14, max = 255),
    rgb(171, 171, 171, max = 255),
    rgb(89, 89, 89, max = 255),
    rgb(95, 158, 209, max = 255),
    rgb(200, 82, 0, max = 255),
    rgb(137, 137, 137, max = 255),
    rgb(162, 200, 236, max = 255),
    rgb(255, 188, 121, max = 255),
    rgb(207, 207, 207, max = 255)
    )
  )

x$parties <- list()

x$parties$BRA <- c(
  PT=rgb(255, 39, 0, max = 255),
  PMDB=rgb(255, 153, 0, max = 255),
  PSDB=rgb(0, 143, 213, max = 255),
  PSB=rgb(213, 94, 0, max = 255),
  PV=rgb(119, 171, 67, max = 255)
  )

x$fte <- c(
  red=rgb(255, 39, 0, max = 255),
  blue=rgb(0, 143, 213, max = 255),
  green=rgb(119, 171, 67, max = 255),
  orange=rgb(230, 159, 0, max = 255)
  )
x$development <- c(
  autumn=rgb(16, 78, 139, max = 255),
  spring=rgb(110, 139, 61, max = 255),
  summer=rgb(154, 50, 205, max = 255),
  winter=rgb(255, 193, 37, max = 255)
)
x$seasons <- c(
  autumn=rgb(16, 78, 139, max = 255),
  spring=rgb(110, 139, 61, max = 255),
  summer=rgb(154, 50, 205, max = 255),
  winter=rgb(255, 193, 37, max = 255)
  )
## return
x
}
NULL


#' @title Color Palettes for Publication (discrete)
#'
#' @description Color palettes for publication-quality graphs. See details.
#' @param palette Palette name.
#'
#' @details The following palettes are available:
#' \itemize{
#' \item {"pub12"}{Default colors of theme_pub}
#' \item {"tableau20"}{Based on software \href{http://www.tableausoftware.com/}{Tableau}}
#' \item {"tableau10"}{Based on software \href{http://www.tableausoftware.com/}{Tableau}}
#' \item {"colorblind"}{Based on software \href{http://www.tableausoftware.com/}{Tableau}}
#'  \item {"tableau10light"}{Based on software \href{http://www.tableausoftware.com/}{Tableau}}
#' }
#' @examples
#' library(scales)
#' show_col(pub_color_pal("pub12")(12))
#' show_col(pub_color_pal("tableau20")(20))
#' show_col(pub_color_pal("tableau10")(10))
#' show_col(pub_color_pal("colorblind")(10))
#' show_col(pub_color_pal("tableau10light")(10))
#' @export
#'
`pub_color_pal` <- function(palette = "pub12") {
  pal.list <- themes_data$pub$colors
  if (!palette %in% c(names(pal.list), "pub12", "tableau10", "tableau20", "tableau10", "tableau10light", "colorblind")) {
    stop(sprintf("%s is not a valid palette name", palette))
  }
  if (palette == "pub12") {
    types <- pal.list[["pub12"]][seq(1, 9, by = 1)]
  } else if (palette == "tableau10") {
    types <- pal.list[["tableau20"]][seq(1, 20, by = 2)]
  } else if (palette == "tableau10light") {
    types <- pal.list[["tableau20"]][seq(2, 20, by = 2)]
  } else if (palette == "colorblind") {
    types <- pal.list[["colorblind"]][seq(1, 10, by = 1)]
  } else {
    types <- pal.list[[palette]]
  }
  function(n) {
    unname(types)[seq_len(n)]
  }
}
NULL


#' @title Publication color scales.
#'
#' @description See \code{\link{pub_color_pal}} for details.
#'
#' @inheritParams ggplot2::scale_colour_hue
#' @inheritParams pub_color_pal
#' @family colour publication
#' @rdname scale_color_pub
#' @export
#' @seealso \code{\link{pub_color_pal}} for references.
#'
scale_colour_pub <- function(palette = "tableau20", ...) {
  discrete_scale("colour", "pub", pub_color_pal(palette), ...)
}
NULL

#' @export
#' @rdname scale_color_pub
scale_fill_pub <- function(palette = "tableau20", ...) {
  discrete_scale("fill", "pub", pub_color_pal(palette), ...)
}
NULL



#' @export
#' @rdname scale_color_pub
scale_color_pub <- scale_colour_pub




#' Extended fivethirtyeight.com color palette
#'
#' The standard fivethirtyeight.com palette for line plots is blue, red, green.
#'  I add an orange ton.
#'
#' @family colour fte
#' @examples
#' library(scales)
#' show_col(fte_color_pal()(4))
#' @export
fte_color_pal <- function() {
  function(n) {
    colors <- themes_data$fte[c("blue", "red", "green", "orange")]
    unname(colors[seq_len(n)])
  }
}
NULL



#' fivethirtyeight.com color scales
#'
#' Color scales using the colors in the fivethirtyeight graphics.
#'
#' @inheritParams ggplot2::scale_colour_hue
#' @family colour fte
#' @rdname scale_fte
#' @seealso \code{\link{theme_538}} for examples.
#' @export
scale_colour_fte <- function(...) {
 discrete_scale("colour", "538", fte_color_pal(), ...)
}
NULL



#' @rdname scale_fte
#' @export
scale_color_fte <- scale_colour_fte

#' @rdname scale_fte
#' @export
scale_fill_fte <- function(...) {
  discrete_scale("fill", "538", fte_color_pal(), ...)
}
NULL




#' @title Color Palettes for Political Organizations (discrete)
#'
#' @description Color palettes for political organizations.
#'
#' @param palette Palette name.
#' @family colour parties
#' @examples
#' library(scales)
#' show_col(parties_color_pal()(10))
#' @export
#'
`parties_color_pal` <- function(palette = "BRA") {
  pal.list <- themes_data$parties
  if (!palette %in% c(names(pal.list), "BRA", "ARG", "CAN", "USA")) {
    stop(sprintf("%s is not a valid palette name", palette))
  }
  if (palette == "BRA") {
    types <- pal.list[["BRA"]][seq(1, 30, by = 1)]
  } else if (palette == "ARG") {
    types <- pal.list[["ARG"]][seq(1, 20, by = 1)]
  } else {
    types <- pal.list[[palette]]
  }
  function(n) {
    unname(types)[seq_len(n)]
  }
}
NULL


#' @title Political Parties Color Scales
#'
#' @description Scale color for political parties.
#'
#' @inheritParams ggplot2::scale_colour_hue
#' @inheritParams parties_color_pal
#' @family colour parties
#' @rdname scale_color_parties
#' @export
#' @seealso \code{\link{parties_color_pal}} for references.
#'
scale_colour_parties <- function(palette = "BRA", ...) {
  discrete_scale("colour", "parties", parties_color_pal(palette), ...)
}
NULL


#' @export
#' @rdname scale_color_parties
scale_fill_parties <- function(palette = "BRA", ...) {
  discrete_scale("fill", "parties", parties_color_pal(palette), ...)
}
NULL



#' @export
#' @rdname scale_color_parties
scale_color_parties <- scale_colour_parties


# http://color.adobe.com/
# http://www.farb-tabelle.de/en/rgb2hex.htm?
# "#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71" "FF7E0D"
