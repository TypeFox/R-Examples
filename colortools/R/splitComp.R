#' @title Split Complementary Color Scheme
#' 
#' @description
#' The split-complementary color scheme is a variation of the complementary
#' color scheme. It uses the two colors adjacent to its complement.
#' 
#' @details
#' This color scheme has the same strong visual contrast as the complementary
#' scheme, but it is supposed to have less tension.
#' 
#' @param color an R color name of a color in hexadecimal notation
#' @param plot logical value indicating whether to plot a color wheel with the
#' generated scheme
#' @param bg background color of the plot. Used only when \code{plot=TRUE}
#' @param labcol color for the labels (i.e. names of the colors). Used only when
#' \code{plot=TRUE}
#' @param cex numeric value indicating the character expansion of the labels
#' @param title logical value indicating whether to display a title in the plot.
#' Unsed only when \code{plot=TRUE}
#' @return A character vector with the given color and the split-complementary
#' colors in hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{complementary}}
#' @export
#' @examples
#' # split-complementary colors of 'tomato' (no plot)
#' splitComp("tomato", plot = FALSE)
#' 
#' # split-complementary colors of 'tomato' (with dark gray background)
#' splitComp("tomato", bg = "gray40")
#'
splitComp <-
function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE)
{	
	tmp_cols = setColors(color, 12)
	split_colors <- tmp_cols[c(1,6,8)]

	# plot
	if (plot)
	{
		# labels color
		if (is.null(labcol)) 
		{
			lab_col = rep("", 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 6, 8)] <- "black"
				lab_col[c(2:5,7,9:12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 6, 8)] <- "white"
				lab_col[c(2:5,7,9:12)] <- col2HSV(bg)
			}
		} else {
			lab_col = rep(labcol, 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 6, 8)] <- labcol
				lab_col[c(2:5,7,9:12)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 6, 8)] <- labcol
				lab_col[c(2:5,7,9:12)] <- col2HSV(bg)
			}
		}	
		# hide non-adjacent colors
		tmp_cols[c(2:5,7,9:12)] <- paste(substr(tmp_cols[c(2:5,7,9:12)],1,7), "0D", sep="")
		pizza(tmp_cols, labcol=lab_col, bg=bg, cex=cex)
		# title
		if (title)
			title(paste("Split complementary colors of: ", tmp_cols[1]), 
				col.main=lab_col[1], cex.main=0.8)
	}
	# result
	split_colors
}
