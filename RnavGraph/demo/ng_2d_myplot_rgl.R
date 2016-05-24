require(RnavGraph) || stop("RnavGraph library not available")
require(rgl) || stop("rgl library not available")
require(MASS) || stop("MASS library not available")

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
				den <- kde2d(x,y)
				
				persp3d(den$x,den$y,den$z, col = "steelblue")  
			}
			
			
			viz1 <- ng_2d_myplot(ng.iris,ng.lg,fnName = "myPlot" , device = "rgl", scaled=TRUE)
			
			nav <- navGraph(ng.iris,ng.lg, viz1)
		})



cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_myplot_rgl.R", package="RnavGraph"),"\n\n\n"))