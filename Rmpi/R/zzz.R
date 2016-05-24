### Copyright (C) 2002 Hao Yu 
.onLoad <- function (lib, pkg) {
    library.dynam("Rmpi", pkg, lib)
    if (!TRUE)
	stop("Fail to load Rmpi dynamic library.")
    if (!is.loaded("mpi_initialize"))
	stop("Probably Rmpi has been detached. Please quit R.")

    if(.Call("mpidist",PACKAGE="Rmpi") == 2){
    	if (length(try(system("lamnodes",TRUE,ignore.stderr = TRUE))) == 0){
    		#cat("\n\tLAM/MPI runtime environment is not operating.\n")
    		#cat("\tStarting LAM/MPI runtime environment.\n")
	    system("lamboot -H",ignore.stderr = TRUE)
	}
    }
	
    if(!.Call("mpi_initialize",PACKAGE = "Rmpi"))
	stop("Cannot start MPI_Init(). Exit")
}

.onUnload <- function(libpath){
    if (mpi.comm.size(1)>0)
		mpi.close.Rslaves()
	mpi.finalize()
	library.dynam.unload("Rmpi", libpath)
}
