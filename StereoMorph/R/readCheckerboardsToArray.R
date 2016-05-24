readCheckerboardsToArray <- function(file, nx, ny, col.reverse = FALSE, row.reverse = FALSE, na.omit=FALSE, ...){

	# CHECK THAT FILE IS MATRIX
	if(!is.matrix(file)) stop("'file' must be a matrix of file paths.")

	# CREATE EMPTY ARRAY
	coor <- array(NA, dim=c(nx*ny, 2, nrow(file), ncol(file)))

	# READ LANDMARKS INTO COORDINATE ARRAY
	for(i in 1:nrow(file)){
		for(j in 1:ncol(file)){

			if(!file.exists(file[i, j])) next
		
			read_table <- as.matrix(read.table(file[i, j], ...))
		
			if(nrow(read_table) != nx*ny) stop(paste0("Dimensions of data in file '", file[i, j], "' do not match input dimensions"))

			coor[, , i, j] <- read_table
		}
	}

	# IF SINGLE SET OF COORDINATES IS INPUT, MAKE SINGLE VALUE MATRIX FOR REVERSE ARRAY
	if(length(dim(coor)) == 2) rev_dim <- c(1, 1)
	
	# IF SINGLE SET OF COORDINATES IS INPUT, MAKE 3D REVERSE ARRAY
	if(length(dim(coor)) == 3) rev_dim <- c(1, 1, dim(coor)[3])

	if(length(dim(coor)) == 4){
		rev_dim <- c(1, 1, dim(coor)[3], dim(coor)[4])

		# IF REVERSE ROW/COLUMN INPUT IS A VECTOR WITH THE SAME LENGTH AS NUMBER OF CAMERA VIEWS MAKE INTO MATRIX TO GET PROPER ORDERING IN ARRAY
		if(is.vector(col.reverse) && length(col.reverse) == dim(coor)[4]) col.reverse <- matrix(col.reverse, nrow=dim(coor)[3], ncol=length(col.reverse), byrow=TRUE)
		if(is.vector(row.reverse) && length(row.reverse) == dim(coor)[4]) row.reverse <- matrix(row.reverse, nrow=dim(coor)[3], ncol=length(row.reverse), byrow=TRUE)
	}

	# MAKE REVERSE ARRAYS
	col_reverse <- array(col.reverse, dim=rev_dim)
	row_reverse <- array(row.reverse, dim=rev_dim)

	# CREATE VECTOR OF INDICES WHERE COLUMNS ARE FLIPPED
	col_ridx <- c(t(matrix(rep(nx:1, ny), nrow=ny, ncol=nx, byrow=TRUE) + matrix(seq(0, nx*ny - nx, by=nx), nrow=ny, ncol=nx)))
	row_ridx <- c(t(matrix(1:(nx*ny), nrow=ny, ncol=nx, byrow=T)[ny:1, ]))

	# FLIP ROWS AND COLUMNS OF EACH CHECKERBOARD POINT ARRAY
	if(length(dim(coor)) == 2){
		if(col_reverse[1, 1]) coor <- coor[col_ridx, ]
		if(row_reverse[1, 1]) coor <- coor[row_ridx, ]
		return(coor)
	}

	if(length(dim(coor)) == 3){
		for(i in 1:dim(coor)[3]) if(col_reverse[, , i]) coor[, , i] <- coor[col_ridx, , i]
		for(i in 1:dim(coor)[3]) if(row_reverse[, , i]) coor[, , i] <- coor[row_ridx, , i]
		return(coor)
	}

	if(length(dim(coor)) == 4){
		for(j in 1:dim(coor)[4]) for(i in 1:dim(coor)[3]) if(col_reverse[, , i, j]) coor[, , i, j] <- coor[col_ridx, , i, j]
		for(j in 1:dim(coor)[4]) for(i in 1:dim(coor)[3]) if(row_reverse[, , i, j]) coor[, , i, j] <- coor[row_ridx, , i, j]
		return(coor)
	}
}