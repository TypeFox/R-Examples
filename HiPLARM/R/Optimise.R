OptimiseAll <- function(increment=128,  verbose=FALSE) {

	print("Optimising Norm for Dge")
	OptimisenormDge(increment, verbose)
	
	print("Optimising rcond for Dge")
	OptimisercondDge(increment, verbose)
	
	print("Optimising determinant for Dge")
	OptimisedetDge(increment, verbose)
	
	print("Optimising crossprod for Dge")
	OptimisecrossprodDge(increment, verbose)
	
	print("Optimising crossprod for Dge,Dge")
	OptimisecrossprodDgeDge(increment, verbose)
	
	print("Optimising crossprod for Dge,matrix")
	OptimisecrossprodDgemat(increment, verbose)	
	
	print("Optimising Solve for Dge")
	OptimiseSolveDge(increment, verbose)
	
	print("Optimising Solve for Dge,matrix")
	OptimiseSolveDgemat(increment, verbose)
	
	print("Optimising %*% for Dge,Dge")
	OptimisematmulDgeDge(increment, verbose)
	
	print("Optimising for LU for Dge")
	OptimiseLU(increment, verbose)	
	
	print("Optimising %*% for Dge,matrix")
	OptimisematmulDgemat(increment, verbose)
	
	print("Optimising Chol for Dpo")	
	OptimiseChol(increment, verbose)
	
	print("Optimising rcond for Dpo")
	OptimisercondDpo(increment, verbose)
	
	print("Optimising solve for Dpo")
	OptimiseSolveDpo(increment, verbose)
	
	print("Optimising Solve for Dpo,Dge")
	OptimiseSolveDpoDge(increment, verbose)
	
	print("Optimising Solve for Dpo,mat")
	OptimiseSolveDpomat(increment, verbose)
	
	print("Optimise %*% for Dtr,Dtr")
	OptimisematmulDtrDtr(increment, verbose)
	
	print("Optimising %*% for Dtr,mat")
	OptimisematmulDtrmat(increment, verbose)
	
	print("Optimising chol2inv for Dtr")
	OptimiseChol2invDtr(increment, verbose)
	
	print("Optimising Solve for Dtr")
	OptimiseSolveDtr(increment, verbose)
	
	print("Optimising Solve for Dtr,mat")
	OptimiseSolveDtrmat(increment, verbose)
	
	hiplarSet("xover_dsyMatrix_matrix_mm", 0)
	hiplarSet("xover_dsyMatrix_norm", 0)
}

checkFile <- function() {

		counter = 0
		
		filename = "default.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))){
			default = scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""), quiet=TRUE)
		} else {
			default = 512
		}


		filename = ".OptnormDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptnormDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""), quiet=TRUE)
		} else {	
			OptnormDge = default
			counter = counter + 1
		}

		filename = ".OptrcondDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptrcondDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""), quiet=TRUE)
		} else {	
			OptrcondDge = default
			counter = counter + 1
		}
		
		filename = ".OptcrossprodDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptcrossprodDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptcrossprodDge = default
			counter = counter + 1
		}
		
		filename = ".OptcrossprodDgeDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptcrossprodDgeDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptcrossprodDgeDge = default
			counter = counter + 1
		}
		
		filename = ".OptcrossprodDgemat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptcrossprodDgemat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptcrossprodDgemat = default
			counter = counter + 1
		}

		filename = ".OptLUDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptLUDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptLUDge = default
			counter = counter + 1
		}

		filename = ".OptdetDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptdetDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptdetDge = default
			counter = counter + 1
		}
	
		filename = ".OptsolveDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDge = default
			counter = counter + 1
		}
		
		filename = ".OptsolveDgemat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDgemat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDgemat = default
			counter = counter + 1
		}

		filename = ".OptmatmulDgemat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptmatmulDgemat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptmatmulDgemat = default
			counter = counter + 1
		}
		
		filename = ".OptrcondDpo.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptrcondDpo <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptrcondDpo = default
			counter = counter + 1
		}

		filename = ".OptsolveDpo.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDpo <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDpo = default
			counter = counter + 1
		}
		
		filename = ".OptsolveDpomat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDpomat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDpomat = default
			counter = counter + 1
		}

		filename = ".OptsolveDpoDge.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDpoDge <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDpoDge = default
			counter = counter + 1
		}
	
		filename = ".OptChol.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptChol <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptChol = default
			counter = counter + 1
		}

		#filename = ".OptnormDtr.dat"
		#if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
		#	OptnormDtr <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		#} else {	
		#	OptnormDtr = default
		#	counter = counter + 1
		#}
	
		filename = ".OptsolveDtr.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDtr <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDtr = default
			counter = counter + 1
		}
		
		filename = ".Optchol2invDtr.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			Optchol2invDtr <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			Optchol2invDtr = default
			counter = counter + 1
		}

		#filename = ".OptsolveDtrDtr.dat"
		#if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/'",filename),collapse=""))) {
		#	OptsolveDtrDtr <- scan(paste(c(find.package("HiPLARM"),"/extdata/'",filename),collapse=""),quiet=TRUE)
		#} else {	
		#	OptsolveDtrDtr = default
		#	counter = counter + 1
		#}
		
		filename = ".OptsolveDtrmat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptsolveDtrmat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptsolveDtrmat = default
			counter = counter + 1
		}
	
		filename = ".OptmatmulDtrmat.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptmatmulDtrmat <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptmatmulDtrmat = default
			counter = counter + 1
		}
	
		filename = ".OptmatmulDtrDtr.dat"
		if(file.exists(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""))) {
			OptmatmulDtrDtr <- scan(paste(c(find.package("HiPLARM"),"/extdata/",filename),collapse=""),quiet=TRUE)
		} else {	
			OptmatmulDtrDtr = default
			counter = counter + 1
		}
		
		.Call("hiplarSet", "xover_dgeMatrix_norm", OptnormDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_rcond", OptrcondDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_crossprod", OptcrossprodDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_dgeMatrix_crossprod", OptcrossprodDgeDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_matrix_crossprod", OptcrossprodDgemat, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_LU", OptLUDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_determinant", OptdetDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_solve", OptsolveDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_matrix_solve", OptsolveDgemat, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dgeMatrix_matrix_mm", OptmatmulDgemat, PACKAGE="HiPLARM")

		.Call("hiplarSet", "xover_dpoMatrix_rcond", OptrcondDpo, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dpoMatrix_solve", OptsolveDpo, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dpoMatrix_matrix_solve", OptsolveDpomat, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dpoMatrix_dgeMatrix_solve", OptsolveDpoDge, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dpoMatrix_chol", OptChol, PACKAGE="HiPLARM")

		.Call("hiplarSet", "xover_dsyMatrix_matrix_mm", 0, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dsyMatrix_norm", 0, PACKAGE="HiPLARM")

		.Call("hiplarSet", "xover_dtrMatrix_solve", OptsolveDtr, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dtrMatrix_chol2inv", Optchol2invDtr, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dtrMatrix_matrix_solve", OptsolveDtrmat, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dtrMatrix_matrix_mm", OptmatmulDtrmat, PACKAGE="HiPLARM")
		.Call("hiplarSet", "xover_dtrMatrix_dtrMatrix_mm", OptmatmulDtrDtr, PACKAGE="HiPLARM")

		print("There are")
		print(counter)
		print("routines not optimised")
		print("For Optimal performance please run the OptimiseAll routine in HiPLARM")

}
