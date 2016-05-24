#' @title Complementary or opposite color
#' 
#' @description
#' Complementary or opposite color scheme is formed by colors that are opposite
#' each other on the color wheel (example: red and green). The high contrast of
#' complementary colors creates a vibrant look that must be managed well so it
#' is not jarring.
#' 
#' @details
#' The complementary color is obtained following a color wheel with 12 colors,
#' each one spaced at 30 degrees from each other.
#' Complementary color schemes are tricky to use in large doses, but work well
#' when you wnat something to stand out. In addition, omplementary colors are
#' really bad for text.
#' 
#' @aliases complementary opposite
#' @param color an R color name or color in hexadecimal notation
#' @param plot logical value indicating whether to plot a color wheel with the
#' generated scheme
#' @param bg background color of the plot. Used only when \code{plot=TRUE}
#' @param labcol color for the labels (i.e. names of the colors). Used only when
#' \code{plot=TRUE}
#' @param cex numeric value indicating the character expansion of the labels
#' @param title logical value indicating whether to display a title in the plot.
#' Used ony when \code{plot=TRUE}
#' @return A character vector with the given color and the complementary color
#' in hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{adjacent}}, \code{\link{splitComp}},
#' \code{\link{triadic}}, \code{\link{tetradic}}, \code{\link{square}}
#' @export
#' @examples
#' # complementary color of 'tomato' with no plot
#' opposite("tomato", plot = FALSE)
#' 
#' # complementary color of 'tomato' with color wheel
#' opposite("tomato", bg = "gray30")
#'
complementary <-
function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE)
{	
	tmp_cols = setColors(color, 12)
	comp_colors <- tmp_cols[c(1, 7)]

	# plot
	if (plot)
	{
		# labels color
		if (is.null(labcol)) 
		{
			lab_col = rep("", 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 7)] <- "black"
				lab_col[c(2:6,8:12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 7)] <- "white"
				lab_col[c(2:6,8:12)] <- col2HSV(bg)
			}
		} else {
			lab_col = rep(labcol, 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 7)] <- labcol
				lab_col[c(2:6,8:12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 7)] <- labcol
				lab_col[c(2:6,8:12)] <- col2HSV(bg)
			}
		}	
		# hide non-adjacent colors
		tmp_cols[c(2:6,8:12)] <- paste(substr(tmp_cols[c(2:6,8:12)],1,7), "0D", sep="")
		pizza(tmp_cols, labcol=lab_col, bg=bg, cex=cex)
		# title
		if (title)
			title(paste("Complementary (opposite) color of: ", tmp_cols[1]), 
				col.main=lab_col[1], cex.main=0.8)
	}
	# result
	comp_colors
}
