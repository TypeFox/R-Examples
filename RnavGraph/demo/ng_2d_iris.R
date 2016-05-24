require(RnavGraph) || stop("RnavGraph library not available")
local({
			## Import the data
			ng.iris <- ng_data(name = "iris", data = iris[,1:4],
					shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
					group = iris$Species,
					labels = substr(iris$Species,1,2))
			
			## get the variable graph node names
			V <- shortnames(ng.iris)
			
			## create the linegraph and its complement
			G <- completegraph(V)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			
			## they are all from the graph
			class(G)
			
			## geberate NG_graph objects
			ng.lg <- ng_graph(name = '3D Transition', graph = LG, layout = 'circle')
			ng.lgnot <- ng_graph(name = '4D Transition', graph = LGnot, layout = 'circle')
			
			## visualization instructions for 2d scatterplots
			viz3dTransition <- ng_2d(ng.iris,ng.lg, glyphs=c("s.L","s.W","p.L","p.W"))
			viz4dTransition  <- ng_2d(ng.iris,ng.lgnot, glyphs=c("s.L","s.W","p.L","p.W"))
			
			## pack them into list
			viz <- list(viz3dTransition, viz4dTransition)	
			graphs <- list(ng.lg, ng.lgnot)
			
			## start navGraph
			nav <- navGraph(data = ng.iris, graph = graphs, viz = viz, settings=list(tk2d=list(linked=FALSE)))
		})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_iris.R", package="RnavGraph"),"\n\n\n"))
