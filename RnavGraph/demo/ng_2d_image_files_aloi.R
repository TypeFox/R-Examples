require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')
require(png) || stop('You need the png package installed!')

local({
			## path to the object images
			imgPath <- system.file("aloi_small", package="RnavGraphImageData")
			aloi_images <- list.files(path=imgPath, full.names=TRUE)
			length(aloi_images)
			aloi_images[1:5]
			
			## sample images
			sel <- sample(1:length(aloi_images),400,replace=FALSE)
			p.aloi_images <- aloi_images[sel]
			
			
			ng.i.objects <- ng_image_files(name="ALOI_objects", path=p.aloi_images)
			ng.i.objects
			
			## import image into R session
			imgData <- t(sapply(p.aloi_images, FUN=function(path){
								x <- readPNG(path)
								r <- sum(x[,,1])
								g <- sum(x[,,2])
								b <- sum(x[,,3])
								return(c(r,g,b))
							}))
			
			
			
			ng.iso.objects <- ng_data(name = "ISO_objects",
					data = data.frame(imgData),
					shortnames = paste('i',1:3, sep = ''))
			
			
			## 3d and 4d transition graph Graphs
			V <- shortnames(ng.iso.objects)
			G <- completegraph(V)
			LG <- linegraph(G)
			LGnot <- complement(LG)
			ng.LG <- ng_graph(name = "3D Transition", graph = LG)
			ng.LGnot <- ng_graph(name = "4D Transition", graph = LGnot)
			
			
			## visualization instructions
			vizObjects1 <- ng_2d(data = ng.iso.objects, graph = ng.LG, images = ng.i.objects)
			vizObjects2 <- ng_2d(data = ng.iso.objects, graph = ng.LGnot, images = ng.i.objects)
			
			
			## start navGraph 
			nav <- navGraph(data = ng.iso.objects, graph = list(ng.LG, ng.LGnot), viz = list(vizObjects1, vizObjects2))
		})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_image_files_aloi.R", package="RnavGraph"),"\n\n\n"))
