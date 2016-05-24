drawLinkage <- function(linkage, method = "svgViewR", file = NULL, animate = TRUE, animate.duration = 1, animate.reverse = FALSE, 
	animate.repeat = -1, path.connect=NULL, connect.joints=TRUE, window.title='Linkage Viewer',
	joint.col.fill="white", joint.col.stroke="black", joint.cex=1.5, joint.lwd=2,
	point.col.fill="black", point.col.stroke="black", point.cex=1, point.lwd=2,
	path.col.fill=NA, path.opacity.fill=1, path.col.stroke="black", path.lwd=1, 
	fade.with.iteration=FALSE, ...){

	#if(is.null(file) && method == "svgViewR") stop("To plot a linkage using the svgViewR method, 'file' must be non-NULL.")
	if(is.null(file)) method <- "plot"

	# FIND FILE TYPE
	if(!is.null(file)){
		extension_found <- FALSE
		if(grepl('[.](jpg$|jpeg$|bmp$|png$|tiff$|eps$)', file, ignore.case=TRUE)){
			method <- "plot"
			extension_found <- TRUE
		}
		if(grepl('[.]html$', file, ignore.case=TRUE)){
			method <- "svgViewR"
			extension_found <- TRUE
		}
		if(grepl('[.][A-Za-z0-9]+$', file, ignore.case=TRUE) && !extension_found) stop("File extension not recognized.")
		if(!extension_found){
			method <- "svgViewR"
			file <- paste0(file, '.html')
		}

		# CHECK FOR FILE WRITING PERMISSION
		tt <- tryCatch(write("", file),error=function(e) e, warning=function(w) w)
		if(!is.null(tt)){
			if(grepl('permission denied', tt, ignore.case=TRUE)){
				stop("It looks like your current working directory is not granting R permission to write files. Please change your current working directory to one with file writing permissions (e.g. Desktop).")
			}
		}
	}

	# COPY POINTS
	points <- linkage$points

	# COPY JOINTS
	joints <- linkage$joint.coor

	# IF MATRIX, CONVERT TO ARRAY
	if(length(dim(joints)) == 2) joints <- array(joints, dim=c(dim(joints)[1], dim(joints)[2], 1))

	# SET PATHS CONNECTING POINTS
	path_list <- NULL
	if(!is.null(points) && !is.null(path.connect)){

		# CREATE PATH LIST TO CONNECT JOINTS
		path_list <- vector("list", length(path.connect))

		if(!is.numeric(path.connect[[1]]) && !is.null(dimnames(points)[[1]])){
			for(i in 1:length(path.connect)){
				for(j in 1:length(path.connect[[i]])){

					# PATH PATH LABELS TO ROWNAMES OF POINT ARRAY
					grepl_match <- grepl(paste0('^', path.connect[[i]][j], '$'), dimnames(points)[[1]])
		
					# ADD TO PATH LIST
					if(sum(grepl_match) > 0) path_list[[i]] <- c(path_list[[i]], which(grepl_match))
				}
			}
		}else{
			for(i in 1:length(path.connect)) path_list[[i]] <- path.connect[[i]]
		}
	}

	if(method == "plot"){

		# GET MATRIX OF ALL POINTS
		all_points <- apply(joints, 2, rbind)
		if(!is.null(points)) all_points <- rbind(all_points, apply(points, 2, rbind))

		# FIND CENTER OF ALL POINTS
		vcenter <- colMeans(all_points, na.rm=TRUE)
		
		# CENTER ALL COORDINATES ABOUT THE ORIGIN
		joints <- joints - array(matrix(vcenter, nrow=dim(joints)[1], ncol=dim(joints)[2], byrow=TRUE), dim=dim(joints))
		if(!is.null(points)) points <- points - array(matrix(vcenter, nrow=dim(points)[1], ncol=dim(points)[2], byrow=TRUE), dim=dim(points))
		all_points <- all_points - matrix(vcenter, nrow=nrow(all_points), ncol=ncol(all_points), byrow=TRUE)
		
		# PROJECT 3D TO 2D
		joints2d <- proj3DTo2D(joints)
		if(!is.null(points)) points2d <- proj3DTo2D(points)
		all_points2d <- proj3DTo2D(all_points)

		# SET GRAPHICAL PARAMETERS
		gp_joints <- c("joint.col.fill", "joint.col.stroke", "joint.cex", "joint.lwd")
		gp_points <- c("point.col.fill", "point.col.stroke", "point.cex", "point.lwd")
		gp_paths <- c("path.col.fill", "path.opacity.fill", "path.col.stroke", "path.lwd")

		# CONVERT GRAPHICAL PARAMETERS TO VECTORS WITH SAME NUMBER OF ELEMENTS OF FIRST X DIMENSION
		for(gpar in gp_joints) if(length(get(gpar)) == 1) assign(gpar, rep(get(gpar), dim(joints)[1]))
		if(!is.null(points)) for(gpar in gp_points) if(length(get(gpar)) == 1) assign(gpar, rep(get(gpar), dim(points)[1]))
		if(!is.null(path.connect)) for(gpar in gp_paths) if(length(get(gpar)) == 1) assign(gpar, rep(get(gpar), length(path.connect)))

		# INITIATE IMAGE FILE
		if(!is.null(file)){
			if(grepl('.bmp$', file, ignore.case=TRUE)) bmp(file, ...)
			if(grepl('.png$', file, ignore.case=TRUE)) png(file, ...)
			if(grepl('.jpg$|.jpeg$', file, ignore.case=TRUE)) jpeg(file, ...)
			if(grepl('.tiff$', file, ignore.case=TRUE)) tiff(file, ...)
			if(grepl('.eps$', file, ignore.case=TRUE)) cairo_ps(file, ...)
		}

		# CREATE PLOT DEVICE THAT FITS ALL JOINT COORDINATES
		par(mar=c(0,0,0,0))
		plot(all_points2d, asp=1, type='n', bty='n', yaxt='n', xaxt='n', xlab='', ylab='')

		for(itr in 1:dim(joints)[3]){

			# PLOT JOINTS
			points(joints2d[, , itr], pch=20, col=joint.col.stroke, cex=joint.cex)

			# CONNECT JOINTS WITH LINES
			if(connect.joints){
				for(i in 1:nrow(linkage$joint.links)){

					# SKIP LINES AMONG JOINTS CONNECTED TO GROUND
					if(linkage$joint.links[i, 'Link.idx'] == 0) next

					# ADD SEGMENTS
					lines(joints2d[c(linkage$joint.links[i, 'Joint1'], linkage$joint.links[i, 'Joint2']), 1, itr], 
						joints2d[c(linkage$joint.links[i, 'Joint1'], linkage$joint.links[i, 'Joint2']), 2, itr])
				}
			}

			# PLOT POINTS
			if(!is.null(points)) points(points2d[, , itr], pch=20, 
				col=point.col.stroke, cex=point.cex)

			# CONNECT POINTS WITH PATHS
			if(!is.null(path_list)){
				for(i in 1:length(path_list)){
				
					path <- path_list[[i]]

					if(is.null(path)) next

					polygon(points2d[path, 1, itr], points2d[path, 2, itr], border=NA, col=path.col.fill[i])

					# DRAW LINES
					if(fade.with.iteration){
						lines(points2d[path, 1, itr], points2d[path, 2, itr], col=rgb(0,0,0, 1 - (itr / (dim(joints)[3]+1))), lwd=path.lwd[i])
					}else{
						lines(points2d[path, 1, itr], points2d[path, 2, itr], col=path.col.stroke[i], lwd=path.lwd[i])
					}
				}
			}
		}
		
		if(!is.null(file)) dev.off()
	}

	if(method == "svgViewR"){
		
		# SET NA TO NONE
		path.col.fill[is.na(path.col.fill)] <- 'none'

		# VECTOR OF DRAWN LINKS
		drawn_links <- c()

		# INITIAL PATH INDEX
		index.add <- 0

		# CREATE NEW SVG FILE
		svgviewr.new(file=file, window.title=window.title, animate.duration=animate.duration, animate.reverse=animate.reverse, animate.repeat)

		# DRAW JOINTS
		if(animate || (length(dim(joints)) > 2 && dim(joints)[3] == 1)){
			svg_points <- svgviewr.points(joints, file=file, col.fill=joint.col.fill, 
				col.stroke=joint.col.stroke, cex=joint.cex, lwd=joint.lwd)
		}
		if(!animate && (length(dim(joints)) > 2 && dim(joints)[3] > 1)){
			for(i in 1:dim(joints)[3]){
				svg_points <- svgviewr.points(joints[, , i], file=file, col.fill=joint.col.fill, 
					col.stroke=joint.col.stroke, cex=joint.cex, lwd=joint.lwd)
			}
		}

		# CONNECT JOINTS WITH PATHS
		if(connect.joints){
			if(animate){
				for(i in 1:nrow(linkage$joint.links)){

					# SKIP LINES AMONG JOINTS CONNECTED TO GROUND
					if(linkage$joint.links[i, 'Link.idx'] == 0) next

					# ADD PATHS
					svgviewr.pathsC(c(linkage$joint.links[i, 'Joint1'], linkage$joint.links[i, 'Joint2']), 
						file=file, index.add=index.add, lwd=path.lwd)
				}
			}else{

				for(itr in 1:dim(joints)[3]){
					
					for(i in 1:nrow(linkage$joint.links)){

						# SKIP LINES AMONG JOINTS CONNECTED TO GROUND
						if(linkage$joint.links[i, 'Link.idx'] == 0) next

						# ADD PATHS
						svgviewr.pathsC(c(linkage$joint.links[i, 'Joint1'], linkage$joint.links[i, 'Joint2']), 
							file=file, index.add=(itr-1)*(nrow(joints)), lwd=path.lwd)
					}
				}
			}
		}

		# ADVANCE POINT COUNT FOR PATH INDEX
		if(animate){
			index.add <- index.add + dim(joints)[1]
		}else{
			index.add <- index.add + (dim(joints)[1]*dim(joints)[3])
		}
		
		if(!is.null(points)){

			# DRAW POINTS
			#svgviewr.points(points, file=file, col.fill=point.col.fill, col.stroke=point.col.stroke, cex=point.cex, lwd=point.lwd)

			if(animate || (length(dim(points)) > 2 && dim(points)[3] == 1)){
				svg_points <- svgviewr.points(points, file=file, col.fill=point.col.fill, 
					col.stroke=point.col.stroke, cex=point.cex, lwd=point.lwd)
			}
			if(!animate && (length(dim(points)) > 2 && dim(points)[3] > 1)){
				for(i in 1:dim(points)[3]){
					svg_points <- svgviewr.points(points[, , i], file=file, col.fill=point.col.fill, 
						col.stroke=point.col.stroke, cex=point.cex, lwd=point.lwd)
				}
			}

			# DRAW PATHS CONNECTING JOINTS
			if(!is.null(path.connect)){

				# CONNECT POINTS WITH PATHS
				path_list_add <- vector("list", length(path.connect))
				if(animate){
					for(j in 1:length(path_list)){
						if(is.null(path_list[[j]])) next
						path_list_add[[j]] <- path_list[[j]] + index.add
					}
					svgviewr.pathsC(path_list_add, file=file, col.fill=path.col.fill, opacity.fill=path.opacity.fill, col.stroke=path.col.stroke, lwd=path.lwd)
				}else{
					for(i in 1:dim(joints)[3]){
						for(j in 1:length(path_list)){
							if(is.null(path_list[[j]])) next
							path_list_add[[j]] <- path_list[[j]] + index.add + dim(points)[1]*(i-1)
						}
						svgviewr.pathsC(path_list_add, file=file, col.fill=path.col.fill, opacity.fill=path.opacity.fill, col.stroke=path.col.stroke, lwd=path.lwd)
					}
				}
			}
		}
	}
	
	1
}
