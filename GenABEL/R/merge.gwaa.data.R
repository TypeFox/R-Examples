"merge.gwaa.data" <-
function(x, y, ... ) {
	if (!is(x,"gwaa.data") || !is(y,"gwaa.data"))
		stop("x and y should have gwaa.data-class")
	snpdata <- try(merge(x@gtdata,y@gtdata, ... ))
	if (is(snpdata,"try-error"))
		stop("error occured in merging gtdata slots of x and y")
	snpdata <- snpdata$data
#	phdata <- merge(x@phdata,y@phdata,by="id",all=T)
	phdata <- merge(x@phdata,y@phdata,all=T)
	if (any(!(snpdata@idnames %in% phdata$id)))
		stop("some ids present in gtdata are missing from phdata")
	if (length(phdata$id) != length(unique(phdata$id)))
		stop("duplicated ids in phdata")
	rownames(phdata) <- phdata$id
	phdata <- phdata[idnames(snpdata),]
	out <- new("gwaa.data",phdata=phdata,gtdata=snpdata)
	out
}
