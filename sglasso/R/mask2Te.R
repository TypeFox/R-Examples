mask2Te <- function(mask){
	e <- mask[upper.tri(mask, diag = FALSE)]
	ge <- setdiff(unique(e), ".")
	ne <- length(ge)
	Te_list <- vector(mode = "list", length = ne)
	Te_ptr_list <- vector(mode = "list", length = ne)
	Te_scn <- vector(mode = "numeric", length = ne + 1)
	Te_ptr_scn <- vector(mode = "numeric", length = ne + 1)
	w <- vector(mode = "numeric", length = ne)
	Te_scn[1] <- 1
	Te_ptr_scn[1] <- 1
	for(i in 1:ne){
		mask_lgCMatrix <- drop0(as(mask == ge[i], "sparseMatrix"))
		Te_list[[i]] <- mask_lgCMatrix@i + 1
		w[i] <- length(mask_lgCMatrix@i)
		Te_scn[i + 1] <- Te_scn[i] + length(Te_list[[i]])
		Te_ptr_list[[i]] <- mask_lgCMatrix@p + 1
		Te_ptr_scn[i + 1] <- Te_ptr_scn[i] + length(Te_ptr_list[[i]])
	}
	Te <- unlist(Te_list)
	Te_ptr <- unlist(Te_ptr_list)
	if(is.character(ge)) nme <- ge
	else nme <- as.character(ge)
	out <- list(Te = Te, Te_scn = Te_scn, Te_ptr = Te_ptr, Te_ptr_scn = Te_ptr_scn, ne = ne, w = w, nme = nme)
	out
}
