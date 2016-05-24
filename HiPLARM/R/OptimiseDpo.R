OptimiseChol <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptChol.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			gc()
			A <- nearPD(Hilbert(size))
			B <- nearPD(Hilbert(size))	
			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(chol(A$mat)))
			else
				timeCPU <- system.time(chol(A$mat))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(chol(B$mat)))
			else
				timeGPU <- system.time(chol(B$mat))

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
	hiplarSet("xover_dpoMatrix_chol",optsize)
}

OptimisercondDpo <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptrcondDpo.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			gc()
			A <- nearPD(Hilbert(size))
			B <- nearPD(Hilbert(size))	
			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(rcond(A$mat)))
			else
				timeCPU <- system.time(rcond(A$mat))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(rcond(B$mat)))
			else
				timeGPU <- system.time(rcond(B$mat))

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
	hiplarSet("xover_dpoMatrix_rcond", optsize)
}

OptimiseSolveDpo <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptsolveDpo.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			gc()
			A <- nearPD(Hilbert(size))
			B <- nearPD(Hilbert(size))	
			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A$mat)))
			else
				timeCPU <- system.time(solve(A$mat))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(B$mat)))
			else
				timeGPU <- system.time(solve(B$mat))

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
	hiplarSet("xover_dpoMatrix_solve", optsize)
}

OptimiseSolveDpoDge <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptsolveDpoDge.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			gc()
			A <- nearPD(Hilbert(size))
			B <- nearPD(Hilbert(size))	
			b <- Matrix(rnorm(size * size), ncol=size)
			a <- Matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A$mat, a)))
			else
				timeCPU <- system.time(solve(A$mat, a))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(B$mat, b)))
			else
				timeGPU <- system.time(solve(B$mat, b))

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
	hiplarSet("xover_dpoMatrix_dgeMatrix_solve", optsize)
}


OptimiseSolveDpomat <- function(increment=128, verbose=FALSE) {

	size <- 256
	time = 0
	complete = 0
	epsilon = 0.0025
	timeout = 0
	optsize = 0
	filename = ".OptsolveDpomat.dat"
	tmp = c(find.package("HiPLARM"),"/extdata/",filename)
	filename = paste(tmp, collapse="")

	if(file.exists(filename)) {

		optsize <- as.numeric(scan(filename))

	} else {
		while(complete == 0 && timeout < 33) {	

			gc()
			A <- nearPD(Hilbert(size))
			B <- nearPD(Hilbert(size))	
			b <- matrix(rnorm(size * size), ncol=size)
			a <- matrix(rnorm(size * size), ncol=size)

			prevtime = time

			hiplarSet("hiplar_library", 1)

			if(verbose)
				timeCPU <- print(system.time(solve(A$mat, a)))
			else	
				timeCPU <- system.time(solve(A$mat, a))

			hiplarSet("hiplar_library", 2)

			if(verbose)
				timeGPU <- print(system.time(solve(B$mat, b)))
			else
				timeGPU <- system.time(solve(B$mat, b))

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
	hiplarSet("xover_dpoMatrix_matrix_solve", optsize)
}
