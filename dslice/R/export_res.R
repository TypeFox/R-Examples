export_res <- function(ds_gsa_obj, file = "", ..., cutoff = 1, decreasing = FALSE, 
                       type = c("name", "size", "DS", "p-value", "FDR", "slice"))
{
	type <- match.arg(type)
	n <- length(ds_gsa_obj$set_name)
	sstr <- vector(length = n, mode = "character")
	snum <- vector(length = n, mode = "numeric")
	for(i in 1:n){
		smat <- ds_gsa_obj$slices[[i]]
		snum[i] <- nrow(smat)
		tmp <- NULL
		for(j in 1:nrow(smat)){
			tmp <- paste(tmp, "| ", smat[j, 1], " + ", smat[j, 2], " = ", smat[j, 3], " ", sep="")
		}
		tmp <- paste(tmp, "|", sep="")
		sstr[i] <- tmp
	}
	dsmat <- cbind(ds_gsa_obj$set_name, ds_gsa_obj$set_size, ds_gsa_obj$DS_value, 
	               ds_gsa_obj$pvalue, ds_gsa_obj$FDR, snum, sstr)
	colnames(dsmat) <- c("gene set", "size", "DS-value", "p-value", "FDR", "|S|", "slices")
	if(type == "name"){
		if(decreasing){
			idx <- order(ds_gsa_obj$set_name, decreasing = TRUE)
		}else{
			idx <- order(ds_gsa_obj$set_name)
		}
	}else if(type == "DS"){
		if(cutoff < 0){
			warning("Threshold is set to be 0.")
			cutoff <- 0
		}
		idx <- which(ds_gsa_obj$DS_value > cutoff)
		if(decreasing){
			idx <- idx[order(ds_gsa_obj$DS_value[idx], decreasing = TRUE)]
		}else{
			idx <- idx[order(ds_gsa_obj$DS_value[idx])]
		}
	}else if(type == "p-value"){
		if((cutoff <= 0) || (cutoff >= 1)){
			warning("Threshold is set to be 1.")
			cutoff <- 1
		}
		idx <- which(ds_gsa_obj$pvalue < cutoff)
		if(decreasing){
			idx <- idx[order(ds_gsa_obj$pvalue[idx], decreasing = TRUE)]
		}else{
			idx <- idx[order(ds_gsa_obj$pvalue[idx])]
		}
	}else if(type == "FDR"){
		if((cutoff <= 0) || (cutoff >= 1)){
			warning("Threshold is set to be 1.")
			cutoff <- 1
		}
		idx <- which(ds_gsa_obj$FDR < cutoff)
		if(decreasing){
			idx <- idx[order(ds_gsa_obj$FDR[idx], decreasing = TRUE)]
		}else{
			dx <- idx[order(ds_gsa_obj$FDR[idx])]
		}
	}else if(type == "slice"){
		if(cutoff < 0){
			warning("Threshold is set to be 1.")
			cutoff <- 1
		}
		idx <- which(snum > cutoff)
		if(decreasing){
			idx <- idx[order(snum[idx], decreasing = TRUE)]
		}else{
			idx <- idx[order(snum[idx])]
		}
	}else if(type == "size"){
		if(decreasing){
			idx <- order(ds_gsa_obj$set_size, decreasing = TRUE)
		}else{
			idx <- order(ds_gsa_obj$set_size)
		}
	}else{
		stop("Invalid type.")
	}
	write.table(dsmat[idx, ], file, ...)
}
