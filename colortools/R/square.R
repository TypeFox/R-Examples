#' @title Square color scheme
#' 
#' @description
#' The square color scheme is similar to the tetradic scheme, but with all four
#' colors spaced around the color circle.
#' 
#' @details
#' The square colors are obtained following a color wheel with 12 colors, each
#' one spaced at 30 degrees from each other. Square color schemes tend to work 
#' best if you let one color be dominant
#' 
#' @param color an R color name or a color in hexadecimal notation
#' @param plot logical value indicating whether to plot a color wheel with the
#' generated scheme
#' @param bg background color of the plot. Used only when \code{plot=TRUE}
#' @param labcol color for the labels (i.e. names of the colors). Used only when
#' \code{plot=TRUE}
#' @param cex numeric value indicating the character expansion of the labels
#' @param title logical value indicating whether to display a title in the plot.
#' Used only when \code{plot=TRUE}
#' @return A character vector with the given color and the square scheme colors
#' in hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{complementary}}, \code{\link{adjacent}},
#' \code{\link{triadic}}, \code{\link{tetradic}}
#' @export
#' @examples
#' # square color scheme for 'tomato'
#' square("tomato")
#'
square <-
function(color, plot=TRUE, bg="white", labcol=NULL, cex = 0.8, title=TRUE)
{	
	tmp_cols = setColors(color, 12)
	sqr_colors <- tmp_cols[c(1,4,7,10)]

	# plot
	if (plot)
	{
		# labels color
		if (is.null(labcol)) 
		{
			lab_col = rep("", 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1,4,7,10)] <- "black"
				lab_col[c(2,3,5,6,8,9,11,12)] <- col2HSV(bg)
			} else {
				lab_col[c(1,4,7,10)] <- "white"
				lab_col[c(2,3,5,6,8,9,11,12)] <- col2HSV(bg)
			}
		} else {
			lab_col = rep(labcol, 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1,4,7,10)] <- labcol
				lab_col[c(2,3,5,6,8,9,11,12)] <- col2HSV(bg)
			} else {
				lab_col[c(1,4,7,10)] <- labcol
				lab_col[c(2,3,5,6,8,9,11,12)] <- col2HSV(bg)
			}
		}	
		# hide non-adjacent colors
		tmp_cols[c(2,3,5,6,8,9,11,12)] <- paste(substr(tmp_cols[c(2,3,5,6,8,9,11,12)],1,7), "0D", sep="")
		pizza(tmp_cols, labcol=lab_col, bg=bg, cex = cex)
		# title
		if (title)
			title(paste("Square color scheme of: ", tmp_cols[1]), 
				col.main=lab_col[1], cex.main=0.8)
	}
	# result
	sqr_colors
}
