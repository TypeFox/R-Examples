"refresh.gwaa.data" <- 
function(data,force=FALSE) {
	if (!is(data,"gwaa.data")) stop("data argument should be of gwaa.data-class")
	err <- try(length(data@gtdata@coding),silent=T)
	if (!is(err,"try-error")) stop("data already appear to be of new format. Use force=T to force")
	data@gtdata@coding <- new("snp.coding",as.raw(rep(1,data@gtdata@nsnps)))
	names(data@gtdata@coding) <- data@gtdata@snpnames
	data@gtdata@strand <- new("snp.strand",as.raw(rep(0,data@gtdata@nsnps)))
	names(data@gtdata@strand) <- data@gtdata@snpnames
	rownames(data@phdata) <- data@phdata$id
	data
}
