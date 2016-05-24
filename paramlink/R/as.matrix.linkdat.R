as.matrix.linkdat = function(x, include.attrs=TRUE, ...) {
	p = do.call(cbind, c(list(FAMID=x$famid, relabel(x$pedigree, x$orig.ids)), x$markerdata))
	if(include.attrs) {
		attr(p, "markerattr") = lapply(x$markerdata, attributes)
      attr(p, "available") = x$available
      attr(p, "model") = x$model
	}
   p
}

restore_linkdat = function(x, attrs=NULL) {
	if(is.null(attrs)) attrs = attributes(x)
	y = linkdat(x[, 1:6, drop=F], model = attrs$model, verbose=FALSE)
	markers = x[, -(1:6), drop=F]
	nMark = ncol(markers)/2
   if(nMark==0) markerdata_list = NULL
	else {
		markerattr = attrs$markerattr
		markerdata_list = lapply(seq_len(nMark), function(k) {
			m = markers[, c(2*k-1,2*k), drop=F]
			attributes(m) = c(markerattr[[k]][-1], list(dim=dim(m)))
			m
		})
		class(markerdata_list) = "markerdata"
	}
	setAvailable(setMarkers(y, markerdata_list), intersect(attrs$available, y$orig.ids))
}

#depreciated
.as.annotated.matrix = function(x) {
	p = cbind(FAMID=x$famid, relabel(x$pedigree, x$orig.ids))
	if(x$nMark>0) {
		p = cbind(p, do.call(cbind, x$markerdata))
		attr(p, "markerattr") = lapply(x$markerdata, attributes)
	}
	attr(p, "available") = x$available
	attr(p, "model") = x$model
	p
}

.restore.linkdat = restore_linkdat