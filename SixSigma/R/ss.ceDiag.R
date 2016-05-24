#' Cause and Effect Diagram
#' 
#' Represents a Cause and Effect Diagram by cause group.
#' 
#' The default value for ss.col is c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD",
#'   "#EEEEEE"), a grayscale style. You can pass any accepted colour string.
#' 
#' @param effect    A short character string that represents the effect we want to analyse.
#' @param causes.gr A vector of characters that represents the causes groups.
#' @param causes    A vector with lists that represents the individual causes for each
#' @param main      Main title for the diagram
#' @param sub       Subtitle for the diagram (recommended the Six Sigma project name)
#' @param ss.col    A vector of colors for a personalized drawing. At least five colors,
#'                   sorted by descendant intensity 
#' @return          A drawing of the causes and effect with "fish-bone" shape
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.\cr
#' 
#' Wikipedia, \url{http://en.wikipedia.org/wiki/Ishikawa_diagram}
#' @note 
#' The cause and effect diagram is also known as "Ishikawa diagram", and
#'   has been widely used in Quality Management. It is one of the Seven
#'   Basic Tools of Quality.
#' @seealso \code{\link{ss.pMap}}
#' @author EL Cano
#' @examples 
#' effect <- "Flight Time"
#' causes.gr <- c("Operator", "Environment", "Tools", "Design", 
#'   "Raw.Material", "Measure.Tool")
#' causes <- vector(mode = "list", length = length(causes.gr))
#' causes[1] <- list(c("operator #1", "operator #2", "operator #3"))
#' causes[2] <- list(c("height", "cleaning"))
#' causes[3] <- list(c("scissors", "tape"))
#' causes[4] <- list(c("rotor.length", "rotor.width2", "paperclip"))
#' causes[5] <- list(c("thickness", "marks"))
#' causes[6] <- list(c("calibrate", "model"))
#' ss.ceDiag(effect, causes.gr, causes, sub = "Paper Helicopter Project")
#' 
#' @export
#' @keywords cause-and-effect
ss.ceDiag <- function(effect, causes.gr, causes, 
		main = "Six Sigma Cause-and-effect Diagram", 
		sub, ss.col = c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")){
    
	n.causes<-length(causes.gr)
	.ss.prepCanvas(main,sub)

#Fish head
	w.head <- unit(1, "strwidth", effect) + unit(4, "mm")
	vp.head <- grid::viewport(x = unit(1, "npc") - w.head,
		height = 1, 
		width = w.head, 
		just = c("left", "center"))
	grid::pushViewport(vp.head)

	x.hshape <- c(0.00, 0.00, 0.50, 0.93, 0.95, 0.93, 0.50)
	y.hshape <- c(0.44, 0.56, 0.55, 0.52, 0.50, 0.48, 0.45)
	grid::grid.xspline(x.hshape, 
		y.hshape, 
		shape = c(0, 0, -1, 1, 0, 1, -1),
		open = FALSE,
	gp = grid::gpar(col = ss.col[2], 
			lwd = 2, fill = ss.col[5]))
	grid::grid.text(effect)
	grid:: popViewport()

#Fish tail
	vp.tail <- grid::viewport(x = 0.01, 
		height = 1, 
		width = 0.05, 
		just = c("left","center"))
	grid::pushViewport(vp.tail)
	grid::grid.xspline(x = c(0, 0, 1),
		y = c(0.44, 0.56, 0.5), 
		shape = c(0, 0, 0), open = FALSE,
		gp = grid::gpar(col = ss.col[2], 
			lwd=2, 
			fill = ss.col[5]))
	grid:: popViewport()
	vp.body <- grid::viewport(x = 0.06,
		height = 1, 
		width = unit(1, "npc") - w.head - unit(0.06, "npc"), 
		just = c("left", "center"))
	grid::pushViewport(vp.body)
	grid:: grid.lines(x = c(0,1), 
		y = c(0.5, 0.5),
		gp = grid::gpar(col = ss.col[2], lwd=4))
	
#body
	#up lines
	m <- (n.causes%/%2) + 1
	pUp <- seq(1/m, 1-(1/m), 1/m)
	for (i in 1:(length(pUp))){
		grid:: grid.lines(x = c(pUp[i] - 0.15,	pUp[i]), 
			y = c(0.8,0.5),
			gp = grid::gpar(col = ss.col[2], lwd=2))
		grid::grid.text(causes.gr[i], 
			x = pUp[i] - 0.15,
			y = 0.81, 
			just = c("center", "bottom"))
		for (j in 1:length(causes[[i]])){
			grid::grid.text(causes[[i]][j], 
				x = unit(pUp[i] - 0.15, "npc") + unit(j, "lines"),
				y = unit(0.80, "npc") - unit(j, "lines"), 
				gp = grid::gpar(fontsize = 8),
				just = c("left", "center"))
			}
		}
	#down lines
	pDown <- pUp + (1/(2*m))
	if (n.causes%%2 != 0) pDown <- c(1/(m*2), pDown)
	k <- length(pDown)
	for (i in (length(pUp)+1):n.causes){
		grid:: grid.lines(x = c((pDown[k] - 0.15), pDown[k]), 
			y = c(0.2, 0.5),
			gp = grid::gpar(col = ss.col[2], lwd=2))
		grid::grid.text(causes.gr[i], 
			x = pDown[k] - 0.15, 
			y = 0.19, 
			just = c("center", "top"))
		for (j in 1:length(causes[[i]])){
			grid::grid.text(causes[[i]][j],
				x = unit((pDown[k] - 0.15), "npc") + unit(j, "lines"),
				y = unit(0.20, "npc") + unit(j, "lines"), 
				gp = grid::gpar(fontsize = 8),
				just=c("left", "center"))
		}
		k <- (k - 1)
	}
}
