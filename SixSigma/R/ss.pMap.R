#' Process Map
#' 
#' This function takes information about the process we want
#'   to represent and draw the Process Map, with its X's, x's, Y's and y's in
#'   each step of the process
#' 
#' The type of the x parameters and y features can be: C(controllable), 
#' N(noise), Cr(Critical), P(Procedure). The default value for ss.col is 
#' c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"), a grayscale style. 
#' You can pass any accepted color string.
#' 
#' 
#' @param steps A vector of characters with the name of the 'n' steps
#' @param inputs.overall A vector of characters with the name of the overall inputs
#' @param outputs.overall A vector of characters with the name of the overall outputs
#' @param input.output A vector of lists with the names of the inputs of the \eqn{i-{th}} step,
#'  that will be the outputs of the \eqn{(i-1)-{th}} step
#' @param x.parameters A vector of lists with a list of the x parameters of the process. The
#'   parameter is a vector with two values: the name and the type (view details)
#' @param y.features   A vector of lists with a list of the y features of the step. The
#'   feature is a vector with two values: the name and the type (view details)
#' @param main The main title for the Process Map
#' @param sub   Subtitle for the diagram (recommended the Six Sigma project name)
#' @param ss.col   A vector of colours for a custom drawing. At least five colours,
#'   sorted by descendant intensity (see details)
#' 
#' @return 
#' A graphic representation of the Map Process.
#' 
#' @references 
#'   \url{http://en.wikipedia.org/wiki/Business_Process_Mapping}\cr
#'  
#'   Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.\cr
#' 
#' @note 
#' The process map is the starting point for a Six Sigma Project, and it is
#'   very important to find out who the x's and y'x are.
#' 
#' @seealso 
#'   \code{\link{ss.ceDiag}}
#' 
#' @author EL Cano
#' 
#' @examples
#' 
#' inputs.overall<-c("operators", "tools", "raw material", "facilities")
#' outputs.overall<-c("helicopter")
#' steps<-c("INSPECTION", "ASSEMBLY", "TEST", "LABELING")
#' #Inputs of process "i" are inputs of process "i+1"
#' input.output<-vector(mode="list",length=length(steps))
#' input.output[1]<-list(c("sheets", "..."))
#' input.output[2]<-list(c("sheets"))
#' input.output[3]<-list(c("helicopter"))
#' input.output[4]<-list(c("helicopter"))
#' 
#' #Parameters of each process
#' x.parameters<-vector(mode="list",length=length(steps))
#' x.parameters[1]<-list(c(list(c("width", "NC")),list(c("operator", "C")),
#' list(c("Measure pattern", "P")), list(c("discard", "P"))))
#' x.parameters[2]<-list(c(list(c("operator", "C")),list(c("cut", "P")),
#' list(c("fix", "P")), list(c("rotor.width", "C")),list(c("rotor.length",
#' "C")), list(c("paperclip", "C")), list(c("tape", "C"))))
#' x.parameters[3]<-list(c(list(c("operator", "C")),list(c("throw", "P")),
#' list(c("discard", "P")), list(c("environment", "N"))))
#' x.parameters[4]<-list(c(list(c("operator", "C")),list(c("label", "P"))))
#' x.parameters
#' 
#' #Features of each process
#' y.features<-vector(mode="list",length=length(steps))
#' y.features[1]<-list(c(list(c("ok", "Cr"))))
#' y.features[2]<-list(c(list(c("weight", "Cr"))))
#' y.features[3]<-list(c(list(c("time", "Cr"))))
#' y.features[4]<-list(c(list(c("label", "Cr"))))
#' y.features
#' 
#' ss.pMap(steps, inputs.overall, outputs.overall,
#'         input.output, x.parameters, y.features, 
#'         sub="Paper Helicopter Project")
#' 
#' @export
#' @keywords process map
ss.pMap <-
function(steps, inputs.overall, outputs.overall,
                  input.output, x.parameters, y.features,
		  main = "Six Sigma Process Map", sub,
		  ss.col = c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")){
	nsteps <- length(steps)
	.ss.prepCanvas(main,sub)
	paintBox <- function(){
		x <- c(rep(0.10, 4), 0.30, 0.70, 
			rep(0.90, 4), 0.7, 0.3)
		y <-c (0.10, 0.37, 0.64, rep(0.91, 4),
			0.64, 0.37, rep(0.1, 3))
		grid::grid.xspline(x, y,
			shape = c(1, 1, 1, 1, 1, 1),
			open = FALSE,
			gp = grid::gpar(fill = ss.col[4], 
				lwd = 3, col = ss.col[2])
			)
	}

#ProcessMap container
	vp.map <- grid::viewport(name = "map", 
			layout.pos.col = 1:3, 
			layout.pos.row = 2,
			layout = grid::grid.layout(3, nsteps, 
					heights = c(0.2, 0.6, 0.2)))	##
	grid::pushViewport(vp.map)

#overall inputs
	vp.inputs <- grid::viewport(layout.pos.col = 1, 
			layout.pos.row = 1, 
			name = "inputs")
	grid::pushViewport(vp.inputs)
	paintBox()
	grid::grid.text("INPUTS\nX")
	grid::grid.move.to(x = 0.5, y = 0.1)
	grid::upViewport()
	vp.inputsText <- grid::viewport(layout.pos.col = 2:nsteps, 
		layout.pos.row = 1, name = "inputst")
	grid::pushViewport(vp.inputsText)
	for (i in 1:length(inputs.overall)){
		grid::grid.text(x = unit(0.5, "cm"),
			y = unit(1, "npc") - unit(i,"lines"),
			paste(inputs.overall[i],"\n"),
			just = c("left","top"), 
			name="inputst")
	}
	grid::upViewport()

#Processes
	for (i in 1:nsteps){
		vp.proc <- grid::viewport(layout.pos.col = i,
			layout.pos.row = 2)
		grid::pushViewport(vp.proc)
		grid::pushViewport(grid::viewport(y = 1, 
				height = 0.5, 
				just = c("center", "top")))
		paintBox()
		grid::grid.lines(x = c(0.1, 0.9), 
			y = c(0.74, 0.74), 
			gp = grid::gpar(lwd = 3, col = ss.col[2]))
		grid::grid.text(steps[i], 
			y = 0.85, 
			just = c("center", "top"))
		grid::grid.text("INPUTS", 
			rot = 90, 
			x = 0.20, 
			y = 0.20,
			just = c("left", "bottom"), gp = grid::gpar(fontsize = 8, col = ss.col[1]))
		for (j in 1:length(input.output[[i]])){
			grid::grid.text(input.output[[i]][j], 
				y = unit(0.7, "npc") - unit(j-1, "lines"),
				just = c("center","top"))
		}
		if (i==1){
			grid::grid.line.to(x = 0.5, y = 0.91,
				arrow = arrow(angle = 30, 
					length = unit(0.15, "inches"),
				ends = "last", 
				type = "open"), 
			gp = grid::gpar(lwd = 6, col = ss.col[1]))
		}
		if (i > 1){
			grid::grid.line.to(x = 0.1, y = 0.5,
				arrow = arrow(angle = 30, length = unit(0.15, "inches"),
				ends = "last", type = "open"), 
				gp=grid::gpar(lwd=6, col=ss.col[1]))
		}
		grid::grid.move.to(x = 0.9, y = 0.5)
		if (i == nsteps){
			grid::grid.move.to(x = 0.7, y = 0.1)
		}
		grid::upViewport()
		grid::pushViewport(grid::viewport(y = 0.5, 
				height = 0.5, 
				just = c("center","top")))
		grid::grid.text("Param.(x): ", 
			y = unit(1,"npc") - unit(2,"mm"), 
			gp = grid::gpar(fontsize = 8),
			x = unit(0,"npc") + unit(2,"mm"),
			just = c("left", "top"))
		for (j in 1:length(x.parameters[[i]])){
			grid::grid.text(paste(x.parameters[[i]][[j]][1], x.parameters[[i]][[j]][2]),
				y = unit(1,"npc") - unit(2,"mm") - unit(j-1,"lines"),
				x = unit(1,"strwidth","Param.(x):   "),
				gp = grid::gpar(fontsize = 8), just = c("left", "top"))
		}
		grid::grid.text("Featur.(y): ",
			y = unit(1, "npc") - 
					unit(2, "mm") - 
					unit(length(x.parameters[[i]]), "lines"),
			x = unit(0, "npc") + unit(2, "mm"),
			gp = grid::gpar(fontsize = 8),
			just = c("left", "top"))
		for (j in 1:length(y.features[[i]])){
			grid::grid.text(y.features[[i]][[j]],
				y = unit(1,"npc") - 
						unit(2,"mm") -
						unit(j - 1 + length(x.parameters[[i]]),	"lines"),
				x = unit(1, "strwidth", "Featur.(y):   "),
				gp = grid::gpar(fontsize = 8), 
				just = c("left", "top"))
		}
		grid::upViewport()
		grid::upViewport()
	}

#overalloutputs
	vp.outputs <- grid::viewport(layout.pos.col = nsteps, 
			layout.pos.row = 3, name = "outputs")
	grid::pushViewport(vp.outputs)
	paintBox()
	grid::grid.text("OUTPUTS\nY")
	grid::grid.line.to(x = 0.7, 
		y = 0.91, 
		arrow = arrow(angle = 30, length = unit(0.15, "inches"),
		ends = "last", type = "open"), 
		gp = grid::gpar(lwd = 6, col = ss.col[1]))
	grid::upViewport()
	vp.outputsText <- grid::viewport(layout.pos.col = 1:(nsteps-1), 
		layout.pos.row = 3, 
		name="outputst")
	grid::pushViewport(vp.outputsText)
	for (i in 1:length(outputs.overall)){
		grid::grid.text(x = unit(1, "npc") - unit(0.5, "cm"),
			y = unit(1, "npc") - unit(i, "lines"), 
			paste(outputs.overall[i],"\n"),
			just = c("right", "top"), 
			name="outputst")
	}
	vp.legend<-grid::viewport(x = unit(0.2, "cm"),
		y = unit(0.2, "cm"),
		just = c("left", "bottom"),
		height = unit(1, "npc") - unit(0.4, "cm"),
		width = 0.3)
	grid::pushViewport(vp.legend)
	grid::grid.rect(gp = grid::gpar(fill = ss.col[3]))
	grid::grid.text("LEGEND\n\t(C)ontrollable\n\t(Cr)itical\n\t(N)oise\n\t(P)rocedure",
		y = unit(1, "npc") - unit(0.2, "cm") , 
		x = unit(0, "npc") + unit(0.2,"cm"),
		just = c("left", "top"),
		gp = grid::gpar(fontsize = 8))
	grid::upViewport()
	grid::upViewport()
}
