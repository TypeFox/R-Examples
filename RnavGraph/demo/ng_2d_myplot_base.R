require(RnavGraph) || stop("RnavGraph library not available")

local({
			ng.iris <- ng_data(name = "iris", data = iris[,1:4],
					shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
					group = iris$Species,
					labels = substr(iris$Species,1,2))
			
			V <- shortnames(ng.iris)
			G <- completegraph(V)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			ng.lg <- ng_graph(name = '3D Transition', graph = LG, layout = 'circle')
			ng.lgnot <- ng_graph(name = '4D Transition', graph = LGnot, layout = 'circle')
			
			
			## Note the plotting functions have to be saved in the .GlobalEnv environment
			## Note, currently order arguments does not work correctly
			myPlot <<- function(x,y,group,labels,order) {
				plot(x,y,col = group, pch = 19)
			}
			
			viz1 <- ng_2d_myplot(ng.iris,ng.lg,fnName = "myPlot" , device = "base", scaled=TRUE)
			
			options(device.ask.default = FALSE) ## omit hit return to see next plot message
			nav <- navGraph(ng.iris,ng.lg, viz1)
			dev.control(displaylist = "inhibit") ## don't svae previous plotting instructions
			
		})



cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_myplot_base.R", package="RnavGraph"),"\n\n\n"))
