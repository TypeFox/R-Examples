require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')

local({
			data(binaryalphadigits)  ## load data

			## create NG_image object
			ng.i.alphadigits <- ng_image_array_gray('BinaryAlphadigits',binaryalphadigits*255,16,20, invert = TRUE)
			
			ng.i.alphadigits  ## preview images
			
			## group label of images
			group <- rep(c(0:9,LETTERS),39)
			
			## isomap, to save time, the isomap coordinates are provided by the RnavGraphImageData package
#			require(vegan) || stop('You need the vegan package installed!')
#			dis <- vegdist(binaryalphadigits)
#			ordalphadigits <- isomap(dis, k=6)
			data(ordalphadigits)
			
			## Data object
			ng.d.alphadigits <- ng_data('BinaryAlphadigits',
					as.data.frame(ordalphadigits$points),
					shortnames = paste('iso',1:dim(ordalphadigits$points)[2], sep = ''),
					group = group)
			
			## Graph
			V <- shortnames(ng.d.alphadigits) ## Node Names
			
			G <- completegraph(V)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			
			## convert these graphs to NG_graph objects
			ng.lg <- ng_graph('3D Transition',LG, layout = 'circle')
			ng.lgnot <- ng_graph('4D Transition',LGnot, layout = 'circle')
			
			## visualize
			viz1 <- ng_2d(graph = ng.lg,
					data = ng.d.alphadigits,
					images = ng.i.alphadigits)
			viz2 <- ng_2d(graph = ng.lgnot,
					data = ng.d.alphadigits,
					images = ng.i.alphadigits)
			
			## Start NavGraph
			nav <- navGraph(data = ng.d.alphadigits,
					graph = list(ng.lg, ng.lgnot),
					viz = list(viz1,viz2))
		})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_images_alpha_letter.R", package="RnavGraph"),"\n\n\n"))


