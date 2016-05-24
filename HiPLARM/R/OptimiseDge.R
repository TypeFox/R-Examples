OptimiseLU <- function(increment=128,verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptLUDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")


	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	
	
				if (size < 1024) {
					timeCPU <- 0
					timeGPU <- 0
					for(i in 1:30 ) {
						A <- Matrix(rnorm(size * size), ncol=size)	
						B <- Matrix(rnorm(size * size), ncol=size)	

							prevtime = time

							hiplarSet("hiplar_library", 1)
							    
								timeCPU <- timeCPU + system.time(lu(A))[3]
								
								hiplarSet("hiplar_library", 2)
							  
								timeGPU <- timeGPU + system.time(lu(B))[3]

					}
							if(verbose) {
								timeCPU <- print(timeCPU)
								timeGPU <  print(timeGPU)
							}	
					diff = (as.numeric(timeGPU)/30) - as.numeric(timeCPU)/30
				} else {

					diff = (as.numeric(timeGPU)) - as.numeric(timeCPU)
			}
				if(verbose)
					print(diff)

						if(diff < 0) {
							optsize = size
								complete = 1
						}


			size = size + increment

				timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_LU", optsize)
}

OptimisecrossprodDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptcrossprodDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")


	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(crossprod(A)))
			else
				timeCPU <- system.time(crossprod(A))

			hiplarSet("hiplar_library", 2)
			
			if(verbose)
				timeGPU <- print(system.time(crossprod(B)))
			else
				timeGPU <- system.time(crossprod(B))
			
			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}

	hiplarSet("xover_dgeMatrix_crossprod", optsize)
}

OptimisecrossprodDgeDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptcrossprodDgeDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")


	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(crossprod(A, B)))
			else
				timeCPU <- system.time(crossprod(A, B))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(crossprod(A, B)))
			else
				timeGPU <- system.time(crossprod(A, B))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_dgeMatrix_crossprod", optsize)

}

OptimisecrossprodDgemat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptcrossprodDgemat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")


	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(crossprod(A, B)))
			else
				timeCPU <- system.time(crossprod(A, B))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(crossprod(A, B)))
			else
				timeGPU <- system.time(crossprod(A, B))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)	

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_matrix_crossprod", optsize)

}

OptimisedetDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptdetDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")


	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(determinant(A)))
			else
				timeCPU <- system.time(determinant(A))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(determinant(B)))
			else
				timeGPU <- system.time(determinant(B))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}

	hiplarSet("xover_dgeMatrix_determinant", optsize)
}

OptimisematmulDgeDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptmatmulDgeDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(A %*% B))
			else
				timeCPU <- system.time(A %*% B)

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(A %*% B))
			else
				timeGPU <- system.time(A %*% B)

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}

	hiplarSet("xover_dgeMatrix_matrix_mm", optsize)
}

OptimisematmulDgemat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptmatmulDgemat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(A %*% B))
			else
				timeCPU <- system.time(A %*% B)

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(A %*% B))
			else
				timeGPU <- system.time(A %*% B)

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_matrix_mm",optsize)
}

OptimisenormDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptnormDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(norm(A, "I")))
			else
				timeCPU <- system.time(norm(A, "I"))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(norm(B, "I")))
			else
				timeGPU <- system.time(norm(B, "I"))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_norm",optsize)

}

OptimiseSolveDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptsolveDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)	
			B <- rnorm(size)

			C <- Matrix(rnorm(size * size), ncol=size)	
			D <- rnorm(size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A,B)))
			else
				timeCPU <- system.time(solve(A,B))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(C,D)))
			else
				timeGPU <- system.time(solve(C,D))
			
			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_solve", optsize)
}

OptimiseSolveDgemat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0	
	optsize = 0
	filename = ".OptsolveDgemat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)	
			B <- matrix(rnorm(size * size), ncol = size)

			C <- Matrix(rnorm(size * size), ncol=size)	
			D <- matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A,B)))
			else
				timeCPU <- system.time(solve(A,B))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(C,D)))
			else
				timeGPU <- system.time(solve(C,D))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + 1
		}
		write(optsize, filename)
	}

	hiplarSet("xover_dgeMatrix_matrix_solve", optsize)
}

OptimisercondDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptrcondDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(rcond(A)))
			else
				timeCPU <- system.time(rcond(A))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(rcond(B)))
			else
				timeGPU <- system.time(rcond(B))

			diff = (as.numeric(timeGPU[3])+epsilon) - as.numeric(timeCPU[3])

			if(verbose)
				print(diff)

			if(diff < 0) {
				optsize = size 
				complete = 1
			}

			size = size + increment

			timeout = timeout + increment
		}
		write(optsize, filename)
	}
	hiplarSet("xover_dgeMatrix_rcond", optsize)
}
