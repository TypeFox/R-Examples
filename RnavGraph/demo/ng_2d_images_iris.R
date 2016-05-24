require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')

local({
			ng.iris <- ng_data(name = "iris", data = iris[,1:4],
					shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
					group = as.numeric(iris$Species),
					labels = substr(iris$Species,1,2))
			
			## Create graph object
			V <- shortnames(ng.iris) ## Node Names
			
			## undirected graph with an adjacency matrix
			adjM <- matrix(c(0,1,1,1, 1,0,1,1, 1,1,0,1, 1,1,1,0),ncol=4, byrow=TRUE)
			G <- newgraph(V,adjM, isAdjacency=TRUE)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			
			## import images
			tclRequire('Img')
			tmpfile <- paste(tempfile(),'.jpg',sep = '')
			download.file('http://www.statlab.uni-heidelberg.de/data/iris/media/Iris-setosa-10_1.jpg',tmpfile)
			id_se <- tclvalue(tkimage.create('photo', file = tmpfile))
			download.file('http://www.statlab.uni-heidelberg.de/data/iris/media/Iris-versicolor-21_1.jpg',tmpfile)
			id_ve <- tclvalue(tkimage.create('photo', file = tmpfile))
			download.file('http://www.statlab.uni-heidelberg.de/data/iris/media/Iris-virginica-3_1.jpg',tmpfile)
			id_vi <- tclvalue(tkimage.create('photo', file = tmpfile))

			ng.img <- new('NG_image', name = 'iris', ids = c(id_se,id_ve,id_vi)[as.numeric(iris$Species)])

	
						
			## convert these graphs to NG_graph objects
			ng.lg <- ng_graph('3D Transition',LG, layout = 'circle')
			ng.lgnot <- ng_graph('4D Transition',LGnot, layout = 'circle')
			
			
			## Visualization
			viz <- list(ng_2d(ng.iris,ng.lg),ng_2d(ng.iris,ng.lg,images = ng.img),ng_2d(ng.iris,ng.lgnot,images=ng.img))
			
			graph <- list(LG = ng.lg, LGnot = ng.lgnot)
			
			## Run navGraph
			nav <- navGraph(ng.iris,graph,viz)
			
		})
cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_images_iris.R", package="RnavGraph"),"\n\n\n"))
