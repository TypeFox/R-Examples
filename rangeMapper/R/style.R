


brewer.pal.get  <- function(palette = NULL) {
	pal = brewer.pal.info
	pal = pal[!pal$category == "qual",]
	bp = lapply(split(pal, row.names(pal)), FUN = function(x) brewer.pal(x$maxcolors, row.names(x)))
	if(!is.null(palette) && palette%in%names(bp) ) bp = bp[palette][[1]]
	bp
	}

#' ggplot theme
#'
#' A ggplot theme based on \code{\link[ggplot2]{theme_bw}}
#'
#' @param base_size   base_size
#' @param base_family base_family
#' @export
#'
theme_rangemap <- function (base_size = 12, base_family = "")  {

    theme_bw(base_size = base_size, base_family = base_family) +
    theme(
		  strip.background     = element_rect(fill   = "grey80",  colour = "grey50", size = 0.2),

		  axis.text            = element_blank(),
		  axis.ticks           = element_blank() ,

		  legend.background    = element_blank(),
		  legend.justification =c(0,0),
		  legend.position      =c(0,0),

		  panel.background     = element_rect(fill   = "white", colour = NA),
		  panel.grid.major     = element_blank(),
		  panel.grid.minor     = element_blank(),
		  panel.border         = element_blank(),

		  plot.title           = element_text(size=  10),
		  plot.background      = element_blank(),
		  plot.margin          = grid::unit(c(1,1,1,1), "mm")

		)
	}

#'
#' A few color palettes
#'
#' @param set   set type (currently one set only)
#' @export
#'
palette_rangemap <- function (set = 'set1')  {

	set1 = c('#5D4893','#38847E','#D7C25B','#D78C5B')
	set1

	}
