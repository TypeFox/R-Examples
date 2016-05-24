#' @title Adjacent or analogous colors
#' 
#' @description
#' Adjacent color schemes use colors that are next to each other on the color
#' wheel. These colors usually match well and create comfortable designs.
#' 
#' @details
#' The analogous colors are obtained following a color wheel with 12 colors,
#' each one spaced at 30 degrees from each other.
#' 
#' @aliases adjacent analogous
#' @param color an R color name or a color in hexadecimal notation
#' @param plot logical value indicating whether to plot a color wheel with the
#' generated scheme
#' @param bg background color of the plot. Used only when \code{plot=TRUE}
#' @param labcol color for the labels (i.e. names of the colors). Used only when
#' \code{plot=TRUE}
#' @param cex numeric value indicating the character expansion of the labels
#' @param title logical value indicating whether to display a title in the plot.
#' Used only when \code{plot=TRUE}
#' @return A character vector with the given color and the analogous colors in
#' hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{complementary}}, \code{\link{splitComp}},
#' \code{\link{triadic}}, \code{\link{tetradic}}, \code{\link{square}}
#' @export
#' @examples
#' # analogous colors of 'red'
#' adjacent("red", plot = FALSE)
#' 
#' # analogous colors of 'tomato' with default color wheel
#' analogous("tomato")
#' 
#' # analogous colors of '#606FEF' with darker background
#' adjacent("#606FEF", bg = "gray20")
#'
adjacent <-
function(color, plot=TRUE, bg="white", labcol=NULL, cex=0.8, title=TRUE)
{	
	tmp_cols = setColors(color, 12)
	adja_colors <- tmp_cols[c(1,2,12)]

	# plot
	if (plot)
	{
		# labels color
		if (is.null(labcol)) 
		{
			lab_col = rep("", 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 2, 12)] <- "black"
				lab_col[c(3:11)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 2, 12)] <- "white"
				lab_col[c(3:11)] <- col2HSV(bg)
			}
		} else {
			lab_col = rep(labcol, 12)
			if (mean(col2rgb(bg)) > 127)
			{
				lab_col[c(1, 2, 12)] <- labcol
				lab_col[c(3:11)] <- col2HSV(bg)
			} else {
				lab_col[c(1, 2, 12)] <- labcol
				lab_col[c(3:11)] <- col2HSV(bg)
			}
		}	
		# hide non-adjacent colors
		tmp_cols[c(3:11)] <- paste(substr(tmp_cols[c(3:11)],1,7), "0D", sep="")
		pizza(tmp_cols, labcol=lab_col, bg=bg, cex=cex)
		# title
		if (title)
			title(paste("Adjacent (analogous) colors of: ", tmp_cols[1]), 
				col.main=lab_col[1], cex.main=0.8)
	}
	# result
	adja_colors
}
