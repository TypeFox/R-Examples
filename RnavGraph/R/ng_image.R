setClass(
		Class = "NG_image",
		representation = representation(
				name = "character",
				ids = "character"
		))


setMethod(f = "show",
		signature = "NG_image",
		definition = function(object){
			
			tt <- tktoplevel()
			tktitle(tt) <- object@name
			fr <- tkframe(tt, borderwidth = 3, relief = 'flat')
			tkpack(fr, side = "top", anchor = "c")
			canvas <- tkcanvas(fr,width = 300, height = 300) #, bg = "lightblue")
			label <- tklabel(fr, text = "Hello")
			n <- length(object@ids)
			slider <- tkscale(fr, from=1, to=n,
					showvalue=T, resolution=1, orient="horizontal", length = 300)
			tkpack(canvas, label, slider, side = "top", anchor = 'c')
			
			
			prev <- tkimage.create('photo')
			
			width <- .tcl2num(tcl('image','width', object@ids[1]))
			height <- .tcl2num(tcl('image','height', object@ids[1]))
			tkconfigure(label, text = paste(width,'x',height, sep = ''))
			maxdim <- max(c(width,height))
			
			if(maxdim > 280) {
				mindim <- min(c(width,height))
				diag <- sqrt(7840*(1+(mindim/maxdim)^2))
				tcl('image_scale', object@ids[1], diag, prev)
			} else {
				tcl(prev,'copy', object@ids[1])
			}
			tkcreate(canvas, 'image', 150, 150, image = prev, tag = 'image')
			
			env <- environment()
			
			showImg <- function(i) {
				tkdelete(canvas,'all')						
				tcl('image', 'delete', prev)
				env$prev <- tkimage.create('photo')
				
				width <- .tcl2num(tcl('image','width', object@ids[i]))
				height <- .tcl2num(tcl('image','height', object@ids[i]))
				tkconfigure(label, text = paste(width,'x',height, sep = ''))
				maxdim <- max(c(width,height))
				
				if(maxdim > 280) {
					mindim <- min(c(width,height))
					diag <- 280*sqrt(1+(mindim/maxdim)^2)
					tcl('image_scale', object@ids[i], diag, prev)
				} else {
					tcl(prev,'copy', object@ids[i])
				}
				tkcreate(canvas, 'image', 150, 150, image = prev, tag = 'image')
			}
			
			tkconfigure(slider, command = function(...){
						i <- as.numeric(...)
						showImg(i)								
					})
		}

)



ng_image_array_gray <- function(name,imageData,width,height, img_in_row = TRUE, invert = FALSE, rotate = 0) {
	
	npix <- width*height
	if(dim(imageData)[ifelse(img_in_row,2,1)] != npix) {
		stop('[ng_image_array_gray]: Data dimension does not match with the width and height argument!')
	}
	
	if(rotate == 0) {
		ii <- 1:npix
		byrow = TRUE
		img_w <- width
		img_h <- height
	}else if(rotate == 90) {
		ii <- 1:npix
		byrow = FALSE
		img_h <- width
		img_w <- height
	}else if(rotate == 180) {
		ii <- npix:1
		byrow = TRUE
		img_w <- width
		img_h <- height
	}else if(rotate == 270) {
		ii <- npix:1
		byrow = FALSE
		img_h <- width
		img_w <- height
	}else {
		stop('[ng_image_array_gray]: Argument rotat must be either 0, 90, 180 or 270!')
	}
	
	## invert
	inv <- ifelse(invert,1,0)
	sign <- ifelse(invert,-1,1)
	
	images <- apply(imageData,
			ifelse(img_in_row,1,2),
			function(img_data) {
				im <- tkimage.create('photo', width = img_w, height = img_h, palette = '256/256/+256')
				tcl(im,'put',
						paste(rbind('{',matrix(grey((inv + sign * img_data/255)[ii]), ncol = img_h, byrow = byrow),'}'),collapse = ' ')
						)
				return(im)
			})
	
	ids <- sapply(images,function(img){tclvalue(img)})
	
	return(new("NG_image", name = name, ids = ids))
}



## TODO: write R functions to import jpgs and pngs if Img tcl extension could not be loaded
ng_image_files <- function(name, paths) {

	unique_paths <- unique(paths)
	unique_ids <- sapply(unique_paths, function(path) {
				tclvalue(tkimage.create('photo', file = path))	
			})
	ii <- match(paths,unique_paths)
	
	return(new("NG_image", name = name, ids = unique_ids[ii]))
}