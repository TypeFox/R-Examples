mask2Tv <- function(mask){
	v <- diag(mask)
	gv <- unique(v)
	nv <- length(gv)
	Tv_rw <- vector(mode = "numeric", length = length(v))
	Tv_pkg <- vector(mode = "numeric", length = length(v))
	Tv_ptr <- vector(mode = "numeric", length = nv + 1)
	Tv_ptr[1] <- 1
	for(i in 1:nv){
		v_ptr <- which(v == gv[i])
		Tv_ptr[i + 1] <- Tv_ptr[i] + length(v_ptr)
		Tv_rw[Tv_ptr[i]:(Tv_ptr[i + 1] - 1)] <- v_ptr
		v_ptr <- v_ptr + v_ptr * (v_ptr - 1)/2
		Tv_pkg[Tv_ptr[i]:(Tv_ptr[i + 1] - 1)] <- v_ptr
	}
	if(is.character(gv)) nmv <- gv
	else nmv <- as.character(gv)
	out <- list(Tv_pkg = Tv_pkg, Tv_rw = Tv_rw,Tv_ptr = Tv_ptr, nv = nv, nmv = nmv)
	out
}