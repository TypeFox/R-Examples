require(RnavGraph) || stop("RnavGraph library not available")
require(RnavGraphImageData) || stop('You need the RnavGraphImageData package installed!')

local({
			data(digits)
			
			## sample some of the data
			sel <- sample(x=1:11000,size = 600)
			p.digits <- digits[,sel]
			group <- rep(c(1:9,0), each = 1100)[sel]
			

ng.i.digits <- ng_image_array_gray('USPS_Handwritten_Digits',digits[,c(1,rep(seq(1100,9900, by = 1100),each = 2)+c(0,1),11000)],16,16,invert = TRUE, img_in_row = FALSE)
ng.i.digits


			## NG_image object
			ng.i.digits <- ng_image_array_gray('USPS_Handwritten_Digits',p.digits,16,16,invert = TRUE, img_in_row = FALSE)
			ng.i.digits
			
			
			
			## isomap
			require(vegan) || stop('You need the vegan package installed!')
			dis <- vegdist(t(p.digits))
			orddigits <- isomap(dis, k=6)
			
			## Data object
			ng.d.digits <- ng_data('USPSHandwrittenDigits',
					as.data.frame(orddigits$points),
					shortnames = paste('iso',1:dim(orddigits$points)[2], sep = ''),
					group = group,
					labels = as.character(group))
			
		
			## start navGraph with certain scagnostics features
			nav <- scagNav(ng.d.digits, scags = "Clumpy", images = ng.i.digits)
			
		})


cat(paste("\n\nThe source code of this demo file is located at:\n",system.file("demo", "ng_2d_images_digits.R", package="RnavGraph"),"\n\n\n"))

