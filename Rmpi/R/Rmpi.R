### Copyright (C) 2002 Hao Yu 
mpi.finalize <- function(){
    if(mpi.is.master())
        print("Exiting Rmpi. Rmpi cannot be used unless relaunching R.")
    .Call("mpi_finalize",PACKAGE = "Rmpi")
}

mpi.exit <- function(){
    if (mpi.is.master())
    	print("Detaching Rmpi. Rmpi cannot be used unless relaunching R.")
    .Call("mpi_finalize",PACKAGE = "Rmpi")
    detach(package:Rmpi)
}

mpi.quit <- function(save="no"){
    .Call("mpi_finalize",PACKAGE = "Rmpi")
    q(save=save,runLast=FALSE)
}

mpi.is.master <- function () 
{
    if (is.loaded("mpi_comm_get_parent"))
	as.logical(.Call("mpi_is_master",PACKAGE = "Rmpi"))
    else {
	if (mpi.comm.size(1)>0)
	    as.logical(mpi.comm.rank(1)==0)
	else
	    as.logical(mpi.comm.rank(0)==0)
    }
}

mpi.any.source <- function(){
    .Call("mpi_any_source",PACKAGE = "Rmpi")
}

mpi.any.tag <- function(){
    .Call("mpi_any_tag",PACKAGE = "Rmpi")
}

mpi.proc.null <- function(){
    .Call("mpi_proc_null",PACKAGE = "Rmpi")
}

string <- function(length){
    if (as.integer(length) < 1)
	stop("need positive length")

    .Call("mkstr",as.integer(length),PACKAGE = "Rmpi")
}

mpi.info.create <- function(info=0){
	.Call("mpi_info_create", as.integer(info),PACKAGE = "Rmpi")
}

mpi.info.set <- function(info=0, key, value){
    .Call("mpi_info_set", as.integer(info), as.character(key), 
	as.character(value),PACKAGE = "Rmpi")
}

mpi.info.get <- function(info=0, key, valuelen){
    .Call("mpi_info_get",as.integer(info), as.character(key), 
	as.integer(valuelen), as.integer(valuelen),PACKAGE = "Rmpi")
}

mpi.info.free <- function(info=0){
	.Call("mpi_info_free", as.integer(info),PACKAGE = "Rmpi")
}

mpi.universe.size <- function(){
	if (!is.loaded("mpi_universe_size")) 
        stop("This function is not supported under MPI 1.2")
	out <-.Call("mpi_universe_size",PACKAGE = "Rmpi")
	if (out==0){
	   # if (exists(".mpi.universe.size"))
		#out<-.mpi.universe.size
	    #else {
			if (.Platform$OS=="windows") {
				require(parallel)
		    	out <- detectCores()
			}
	    #}		
	}
	if (.Call("mpidist",PACKAGE="Rmpi") == 2)
	    out <- out-length(grep("no_schedule",system("lamnodes",TRUE,ignore.stderr=TRUE)))
	if (.Call("mpidist",PACKAGE="Rmpi") == 1 && out == 1){
		if (length(unlist(strsplit(.Platform$pkgType,"mac"))) ==2)
			out <- as.integer(unlist(strsplit(system("sysctl hw.ncpu",TRUE,ignore.stderr=TRUE),":"))[2])
	}
	#if (.Call("mpidist",PACKAGE="Rmpi") == 1 && out > 1)
	#	if (.Platform$OS!="windows")
	#		out <- out-1
	out
}

mpi.get.processor.name <- function(short=TRUE){
    name <- .Call("mpi_get_processor_name",PACKAGE = "Rmpi")
    if (short)
	name <- unlist(strsplit(name, "\\."))[1]
    name
}

mpi.sendrecv <-  function(senddata, sendtype, dest, sendtag, recvdata, 
			recvtype, source, recvtag, 
         		comm = 1, status = 0) 
 {
   .Call("mpi_sendrecv", senddata, as.integer(sendtype), 
	  as.integer(dest), 
          as.integer(sendtag), recvdata, as.integer(recvtype), 
          as.integer(source), as.integer(recvtag), as.integer(comm),
          as.integer(status), PACKAGE="Rmpi")
}

mpi.sendrecv.replace <- function(x, type, dest, sendtag, source, recvtag,  
         comm = 1, status = 0)
 {
   .Call("mpi_sendrecv_replace", x, as.integer(type), as.integer(dest),
          as.integer(sendtag), as.integer(source), as.integer(recvtag), 
          as.integer(comm), as.integer(status), PACKAGE="Rmpi")
}

mpi.cart.create <- function(commold=1, dims, periods, reorder=FALSE, commcart=3) {
        .Call("mpi_cart_create", as.integer(commold), as.integer(dims), 
        as.integer(periods), as.integer(reorder), as.integer(commcart), PACKAGE="Rmpi")
}

mpi.cartdim.get <- function(comm=3) {
        .Call("mpi_cartdim_get",as.integer(comm), PACKAGE="Rmpi")
}

mpi.cart.get <- function(comm=3, maxdims) {

        out <- .Call("mpi_cart_get",as.integer(comm), as.integer(maxdims), PACKAGE="Rmpi")
        dims <- out[1:maxdims]
        periods <- out[(maxdims+1):(maxdims*2)]
        coords <- out[(maxdims*2 + 1):(maxdims*3)]
        list(dims=dims,periods=periods,coords=coords)
}

mpi.cart.rank <- function(comm=3, coords) {
        .Call("mpi_cart_rank",as.integer(comm), as.integer(coords), PACKAGE="Rmpi")
}

mpi.cart.coords <- function(comm=3, rank, maxdims) {
        .Call("mpi_cart_coords",as.integer(comm), as.integer(rank), as.integer(maxdims), 
	PACKAGE="Rmpi")
}

mpi.cart.shift <- function(comm=3, direction, disp){
	.Call("mpi_cart_shift",   as.integer(comm), as.integer(direction-1), 
		as.integer(disp), PACKAGE="Rmpi")
}

mpi.dims.create <- function(nnodes, ndims, dims=integer(ndims)){
	.Call("mpi_dims_create",as.integer(nnodes),as.integer(ndims),as.integer(dims),
	PACKAGE="Rmpi")
}
