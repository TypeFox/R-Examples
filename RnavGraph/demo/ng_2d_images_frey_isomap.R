require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')

local({
			data(frey)
			## isomap data dimensionality reduction
			#require(vegan) || stop("library vegan is needed for this demo.")

			#frey2 <- t(frey)
			dims <- 6 
			#dise <- vegdist(frey2, method="euclidean")
			#ord <- isomap(dise,k = 12, ndim= dims, fragmentedOK = TRUE)
			data(ordfrey)
			
			iso.frey <- data.frame(ordfrey$points)
			
	
			## Images
			## sample a few
			sel <- seq(1,dim(frey)[2],3)
			ng.frey <- ng_image_array_gray('Brendan_Frey',frey[,sel],28,20, img_in_row = FALSE, rotate = 90)
			ng.frey
			
			
			V <- mapply(function(x){paste("i",x,sep="")},1:dims)
			ng.iso.frey <- ng_data(name = "ISO_frey",
						data = iso.frey[sel,],
						shortnames = V)
				
			#nav <- scagNav(ng.iso.frey,c("Clumpy","Outlying","Skewed","Monotonic", "Striated"),
			#				images=ng.frey)
			
			G <- completegraph(V)
			
			LG <- linegraph(G)
			LGnot <- complement(LG)
			
			ng.LG <- ng_graph("3d frey", LG)
			ng.LGnot <- ng_graph("4d frey", LGnot)
			
			viz3d <- ng_2d(ng.iso.frey,ng.LG, images = ng.frey)
			viz4d <- ng_2d(ng.iso.frey,ng.LGnot, images = ng.frey)
			
			nav <- navGraph(ng.iso.frey,list(ng.LG,ng.LGnot),list(viz3d,viz4d))
			
			
		})

cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_images_frey_isomap.R", package="RnavGraph"),"\n\n\n"))
