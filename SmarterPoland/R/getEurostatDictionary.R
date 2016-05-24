getEurostatDictionary <-
	function(dictname) {
  read.table(paste("http://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=dic%2Fen%2F",dictname,".dic",sep=""), sep="\t", header=F, stringsAsFactors=FALSE)
}
