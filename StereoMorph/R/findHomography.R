findHomography <- function(coor.img, coor.obj=NULL, nx=NULL, ny=NULL){

	## 11/11/2015 I checked results checked against opencv code results and they were the same
	
	if(is.null(coor.obj) && is.null(nx) && is.null(ny)) stop(paste0('coor.obj is NULL, either nx or ny must be specified.'))

	if(is.null(coor.obj)){
	
		if(is.null(nx)) nx <- dim(coor.img)[1] / ny
		if(is.null(ny)) ny <- dim(coor.img)[1] / nx

		# CREATE DEFAULT OBJECT COORDINATE MATRIX
		coor.obj <- cbind(
			c(matrix(t(matrix(0:(ny-1), nrow=ny, ncol=nx)), nrow=1, byrow=F)), 
			rep(0:(nx-1), ny), 
			rep(0, nx*ny)
		)

		# FILL ARRAY OF SAME DIMENSIONS AS 2D COORDINATES
		if(length(dim(coor.img)) > 2) coor.obj <- array(coor.obj, dim = c(dim(coor.obj), dim(coor.img)[3]))
	}

	# CONVERT ALL INPUT DIMENSIONS TO 4D
	if(length(dim(coor.img)) == 2){
		idim <- 2
		coor.obj <- array(coor.obj, dim=c(dim(coor.obj), 1, 1), dimnames=list(dimnames(coor.obj)[[1]], dimnames(coor.obj)[[2]], NULL, NULL))
		coor.img <- array(coor.img, dim=c(dim(coor.img), 1, 1), dimnames=list(dimnames(coor.img)[[1]], dimnames(coor.img)[[2]], NULL, NULL))
	} else if(length(dim(coor.img)) == 3){
		idim <- 3
		coor.obj <- array(coor.obj, dim=c(dim(coor.obj), 1), dimnames=list(dimnames(coor.obj)[[1]], dimnames(coor.obj)[[2]], dimnames(coor.obj)[[3]], NULL))
		coor.img <- array(coor.img, dim=c(dim(coor.img), 1), dimnames=list(dimnames(coor.img)[[1]], dimnames(coor.img)[[2]], dimnames(coor.img)[[3]], NULL))
	} else if(length(dim(coor.img)) == 4) idim <- 4

	H <- array(NA, dim=c(3, 3, dim(coor.img)[3:4]), dimnames=list(NULL, NULL, dimnames(coor.img)[[3]], dimnames(coor.img)[[4]]))
	error <- array(NA, dim=c(dim(coor.img)[1],dim(coor.img)[3:4]), dimnames=list(dimnames(coor.img)[[1]], dimnames(coor.img)[[3]], dimnames(coor.img)[[4]]))

	for(k in 1:dim(coor.img)[3]){
		for(m in 1:dim(coor.img)[4]){

			# EMPTY COEFFICIENT MATRIX
			a <- matrix(NA, nrow=2*nrow(coor.img[, , k, m]), ncol=9)
		
			# FILL COEFFICIENT MATRIX
			j <- 1
			for(i in seq(1, nrow(coor.img)*2, by=2)){
				a[i, ] <- c(coor.obj[j, 1, k, m], coor.obj[j, 2, k, m], 1, 0, 0, 0, -coor.img[j, 1, k, m]*coor.obj[j, 1, k, m], -coor.img[j, 1, k, m]*coor.obj[j, 2, k, m], -coor.img[j, 1, k, m])
				a[i+1, ] <- c(0, 0, 0, coor.obj[j, 1, k, m], coor.obj[j, 2, k, m], 1, -coor.img[j, 2, k, m]*coor.obj[j, 1, k, m], -coor.img[j, 2, k, m]*coor.obj[j, 2, k, m], -coor.img[j, 2, k, m])
				j <- j + 1
			}
		
			# SOLVE USING EIGEN ANALYSIS
			eigen_r <- eigen(t(a) %*% a)
		
			# GET EIGENVECTOR CORRESPONDING TO MINIMUM EIGENVALUE
			H[, , k, m] <- eigen_r$vectors[, which.min(eigen_r$values)]
			
			# CONVERT VECTOR TO MATRIX
			H[, , k, m] <- matrix(H[, , k, m], nrow=3, byrow=TRUE)

			# FIND INVERSE
			H[, , k, m] <- solve(H[, , k, m])

			# SCALE HOMOGRAPHY MATRIX SO THAT 3,3 = 1
			H[, , k, m] <- H[, , k, m] * (1 / H[3,3,k,m])
		
			# FIND ERROR IN HOMOGRAPHY (OBJECT INVERSE)
			#coor_obj_inv <- t(H[, , k, m] %*% t(cbind(coor.img[, , k, m], 1)))
			#error[, k, m] <- sqrt(rowSums((coor_obj_inv - cbind(coor.obj[, 1:2, k, m], 1))^2))

			# FIND ERROR IN HOMOGRAPHY (IMAGE INVERSE - BEST BECAUSE IN STANDARD PIXEL COORDINATES)
			coor_img_inv <- t(solve(H[, , k, m]) %*% t(cbind(coor.obj[, 1:2, k, m], 1)))
			coor_img_inv <- coor_img_inv * (1 / coor_img_inv[, 3])
			error[, k, m] <- sqrt(rowSums((coor_img_inv - cbind(coor.img[, 1:2, k, m], 1))^2))
			#print(range(error[, k, m]))
		}
	}

	# RETURN HOMOGRAPHY MATRIX OR ARRAY OF DIMENSIONS CORRESPONDING TO INPUT
	if(idim == 2){
		rlist <- list(
			H = H[, , 1, 1],
			error = error[, 1, 1]
		)
		return(rlist)
	}
	if(idim == 3){
		rlist <- list(
			H = H[, , , 1],
			error = error[, , 1]
		)
		return(rlist)
	}
	rlist <- list(
		H = H,
		error = error
	)
}