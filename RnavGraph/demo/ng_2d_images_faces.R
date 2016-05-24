require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')
require(PairViz) || stop('You need the PairViz package installed!')

local({
			data(faces)
			ng.faces <- ng_image_array_gray('Olivetti Faces',faces,64,64, img_in_row = FALSE)
			ng.faces
			
			group <- rep(LETTERS[1:40], each = 10)
			
			## data
			D <- L2Distance(as.matrix(faces),as.matrix(faces))
			Coord <- cmdscale(D,k=4)
			
			## principal coordinate analysis
			ng.pca <- ng_data('olivetti', data = as.data.frame(Coord), shortnames = paste('pc',1:4,sep = ''), group = group)
			
			
			## graphs
			V <- shortnames(ng.pca)
			G <- completegraph(V)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			
			ng.LG <- ng_graph("3d Transition",LG)
			ng.LGnot <- ng_graph("4d Transition",LGnot)
			
			viz1 <- ng_2d(ng.pca, ng.LG, images = ng.faces, glyphs = eulerian(as(G,"graphNEL")))
			viz2 <- ng_2d(ng.pca, ng.LG, images = ng.faces, glyphs = eulerian(as(G,"graphNEL")))
			
			navGraph(data = ng.pca, graph = ng.LG, viz = list(viz1,viz2), settings = NULL)

		})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_images_faces.R", package="RnavGraph"),"\n\n\n"))
