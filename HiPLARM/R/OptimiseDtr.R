OptimiseChol2invDtr <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	filename = ".Optchol2invDtr.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	optsize = 0

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 45) {	


			A <- Matrix(rnorm(size * size), ncol=size)	
			A <- triu(A)
			B <- rnorm(size)

			C <- Matrix(rnorm(size * size), ncol=size)	
			C <- triu(C)
			D <- rnorm(size)

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
	hiplarSet("xover_dtrMatrix_dtrMatrix_mm", optsize)
}

OptimisematmulDtrDtr <- function(increment=128, verbose=FALSE) {

	  size <- 256
	  time = 0
	  complete = 0
	  epsilon = 0.0025
	  timeout = 0
	  optsize = 0
	  filename = ".OptmatmulDtrDtr.dat"
	  tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	  filename = paste(tmp, collapse="")

	  if(file.exists(filename)) {

			optsize <- as.numeric(scan(filename))

	  } else {
			while(complete == 0 && timeout < 33) {

				  A <- Matrix(rnorm(size * size), ncol=size)
				  B <- Matrix(rnorm(size * size), ncol=size)
				  A <- triu(A)
				  B <- triu(B)

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
	  hiplarSet("xover_dtrMatrix_dtrMatrix_mm", optsize)
}

OptimisematmulDtrmat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptmatmulDtrmat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			A <- Matrix(rnorm(size * size), ncol=size)
			a <- matrix(rnorm(size * size), ncol=size)
			B <- Matrix(rnorm(size * size), ncol=size)
			b <- matrix(rnorm(size * size), ncol=size)
			A <- triu(A)		
			B <- triu(B)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(A %*% a))
			else
				timeCPU <- system.time(A %*% a)

			hiplarSet("hiplar_library", 2)


			if(verbose)
				timeGPU <- print(system.time(B %*% b))
			else
				timeGPU <- system.time(B %*% b)

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
	hiplarSet("xover_dtrMatrix_matrix_mm", optsize)
}

OptimiseSolveDtr <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	filename = ".OptsolveDtr.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	optsize = 0

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 45) {	


			A <- Matrix(rnorm(size * size), ncol=size)	
			A <- triu(A)
			B <- rnorm(size)

			C <- Matrix(rnorm(size * size), ncol=size)	
			C <- triu(C)
			D <- rnorm(size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A)))
			else	
				timeCPU <- system.time(solve(A))


			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(C)))
			else
				timeGPU <- system.time(solve(C))

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
	hiplarSet("xover_dtrMatrix_solve", optsize)
}

OptimiseSolveDtrmat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	filename = "../extdata/.OptsolveDtrmat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	optsize = 0

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))


	} else {
		while(complete == 0 && timeout < 45) {	


			A <- Matrix(rnorm(size * size), ncol=size)	
			A <- triu(A)
			B <- matrix(rnorm(size * size), ncol=size)

			C <- Matrix(rnorm(size * size), ncol=size)	
			C <- triu(C)
			D <- matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A, B)))
			else
				timeCPU <- system.time(solve(A, B))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(C, D)))
			else
				timeGPU <- system.time(solve(C, D))

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
	hiplarSet("xover_dtrMatrix_matrix_solve", optsize)
}
