get.meta.info <-
function(app.id, stat.id){
	meta.info.api.uri <- get.api.url(app.id, 'meta-info', list(), stat.id)
	meta.info <- scan(meta.info.api.uri, list(obj.id="", obj.name="", name="", code="", level="", unit=""), sep="\t")
	return(as.data.frame(meta.info))
}

