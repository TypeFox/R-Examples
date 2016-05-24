#' @title Triadic Color Scheme
#' 
#' @description
#' Triadic color schemes use colors that are evenly spaced around the color
#' wheel.
#' 
#' @details
#' The triadic colors are obtained following a color wheel with 12 colors, each
#' one spaced at 30 degrees from each other.
#' Triadic color schemes tend to be quite vibrant. To use a triadic harmony
#' successfully, the colors should be carefully balanced letting one color
#' dominate and use the others for accent.
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
#' @return A character vector with the given color and the triadic colors in
#' hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{complementary}}, \code{\link{splitComp}},
#' \code{\link{adjacent}}, \code{\link{tetradic}}, \code{\link{square}}
#' @export
#' @examples
#' # triadic colors of 'tomato'
#' triadic("tomato")
#' 
#' # triadic colors of 'tomato' with background color 'gray20'
#' triadic("tomato", bg = "gray20")
#'
triadic <-
function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE)
{	
	tmp_cols = setColors(color, 12)
	triad_colors <- tmp_cols[c(1,5,9)]

	# plot
	if (plot)
	{
		# labels color
		if (is.null(labcol)) 
		{
			lab_col = rep("", 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 5, 9)] <- "black"
				lab_col[c(2,3,4,6,7,8,10,11,12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 5, 9)] <- "white"
				lab_col[c(2,3,4,6,7,8,10,11,12)] <- col2HSV(bg)
			}
		} else {
			lab_col = rep(labcol, 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 5, 9)] <- labcol
				lab_col[c(2,3,4,6,7,8,10,11,12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 5, 9)] <- labcol
				lab_col[c(2,3,4,6,7,8,10,11,12)] <- col2HSV(bg)
			}
		}
		# hide non-adjacent colors
		tmp_cols[c(2,3,4,6,7,8,10,11,12)] <- paste(
			substr(tmp_cols[c(2,3,4,6,7,8,10,11,12)],1,7), "0D", sep="")
		pizza(tmp_cols, labcol=lab_col, bg=bg, cex=cex)
		# title
		if (title)
			title(paste("Triadic color scheme of: ", tmp_cols[1]), 
				col.main=lab_col[1], cex.main=0.8)
	}
	# result
	triad_colors
}
