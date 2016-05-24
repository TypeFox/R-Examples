sglasso.fit <- function(Sv, mask, w, flg, nrho, min_rho, nstep, algorithm, truncate, tol){
	out_mask2Tv <- mask2Tv(mask)
	Tv_pkg <- out_mask2Tv$Tv_pkg
	Tv_rw <- out_mask2Tv$Tv_rw
	Tv_ptr <- out_mask2Tv$Tv_ptr
	nv <- out_mask2Tv$nv
	nmv <- out_mask2Tv$nmv
	out_mask2Te <- mask2Te(mask)
	Te <- out_mask2Te$Te
	Te_scn <- out_mask2Te$Te_scn
	Te_ptr <- out_mask2Te$Te_ptr
	Te_ptr_scn <- out_mask2Te$Te_ptr_scn
	ne <- out_mask2Te$ne
	nme <- out_mask2Te$nme
	nSv <- length(Sv)
	nTv <- Tv_ptr[nv + 1] - 1
	nTe <- Te_scn[ne + 1] - 1
	nTe_ptr <- Te_ptr_scn[ne + 1] - 1
	if(is.null(w)) w <- out_mask2Te$w
	else{
		if(any(w <= 0)) stop("weights must be non-negative values")
		if(length(w) != ne) stop("length of 'w' is not equal to 'ne'")
	}
	if(is.null(flg)) flg <- rep(1, ne)
	else{
		if(length(flg) != ne) stop("length of 'flg' is not equal to 'ne'")
		if(!is.logical(flg)) stop("flg is not a logical vector")
	}
	storage.mode(nSv) <- "integer"
	storage.mode(Sv) <- "double"
	storage.mode(nTv) <- "integer"
	storage.mode(Tv_pkg) <- "integer"
	storage.mode(Tv_rw) <- "integer"
	storage.mode(nv) <- "integer"
	storage.mode(Tv_ptr) <- "integer"
	storage.mode(nTe) <- "integer"
	storage.mode(Te) <- "integer"
	storage.mode(nTe_ptr) <- "integer"
	storage.mode(Te_ptr) <- "integer"
	storage.mode(ne) <- "integer"
	storage.mode(Te_scn) <- "integer"
	storage.mode(Te_ptr_scn) <- "integer"
	storage.mode(nstep) <- "integer"
    storage.mode(truncate) <- "double"
	storage.mode(tol) <- "double"
	rho <- double(nrho)
	storage.mode(nrho) <- "integer"
	storage.mode(min_rho) <- "double"
	grd <- matrix(0, nv + ne, nrho, dimnames = list(c(nmv, nme)))
	storage.mode(grd) <- "double"
	th <- matrix(0, nv + ne, nrho, dimnames = list(c(nmv, nme)))
	storage.mode(th) <- "double"
	storage.mode(w) <- "double"
	storage.mode(flg) <- "integer"
	n <- integer(1)
	conv <- integer(1)
	if(nrho > 1){
		if(algorithm == "ccd"){
			fit = .Fortran("sglasso_ccd_path", nSv = nSv, Sv = Sv, nTv = nTv, Tv_pkg = Tv_pkg, Tv_rw = Tv_rw, nv = nv, 
						   Tv_ptr = Tv_ptr, nTe = nTe, Te = Te, nTe_ptr = nTe_ptr, Te_ptr = Te_ptr,ne = ne, 
						   Te_scn = Te_scn, Te_ptr_scn = Te_ptr_scn, nstep = nstep, trnc = truncate, tol = tol, rho = rho,
						   nrho = nrho, min_rho = min_rho, grd = grd, th = th, w = w, pnl_flg = flg, n = n, conv = conv)
		} else {
			fit = .Fortran("sglasso_ccm_path", nSv = nSv, Sv = Sv, nTv = nTv, Tv_pkg = Tv_pkg, Tv_rw = Tv_rw, nv = nv, 
						   v_ptr = Tv_ptr, nTe = nTe, Te = Te, nTe_ptr = nTe_ptr, Te_ptr = Te_ptr,ne = ne, 
						   Te_scn = Te_scn, Te_ptr_scn = Te_ptr_scn, nstep = nstep, trnc = truncate, tol = tol, rho = rho, 
						   nrho = nrho, min_rho = min_rho, grd = grd, th = th, w = w, pnl_flg = flg, n = n, conv = conv)
		}
	} else {
		if(algorithm == "ccd"){
			fit = .Fortran("sglasso_ccd_single", nSv = nSv, Sv = Sv, nTv = nTv, Tv_pkg = Tv_pkg, Tv_rw = Tv_rw, nv = nv, 
						   Tv_ptr = Tv_ptr, nTe = nTe, Te = Te, nTe_ptr = nTe_ptr, Te_ptr = Te_ptr,ne = ne, 
						   Te_scn = Te_scn, Te_ptr_scn = Te_ptr_scn, nstep = nstep, trnc = truncate, tol = tol, rho = min_rho,
						   grd = grd, th = th, w = w, pnl_flg = flg, n = n, conv = conv)
		} else {
			fit = .Fortran("sglasso_ccm_single", nSv = nSv, Sv = Sv, nTv = nTv, Tv_pkg = Tv_pkg, Tv_rw = Tv_rw, nv = nv, 
						   Tv_ptr = Tv_ptr, nTe = nTe, Te = Te, nTe_ptr = nTe_ptr, Te_ptr = Te_ptr,ne = ne, 
						   Te_scn = Te_scn, Te_ptr_scn = Te_ptr_scn, nstep = nstep, trnc = truncate, tol = tol, rho = min_rho,
						   grd = grd, th = th, w = w, pnl_flg = flg, n = n, conv = conv)			
		}
	}
	if(fit$conv!=0) warning("sglasso does not converge with code ", fit$conv, "\n")
	fit
}
