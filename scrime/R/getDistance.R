`getDistance` <-
function(data){		
	tmp.le<-length(data)
	tmp.samp<-sample(data,min(tmp.le,1000))
	tmp.le<-length(unique(tmp.samp))
	distance<-ifelse(tmp.le>10,"euclidean","smc")
	warning("Since distance has not been specified, ",distance," is used as distance measure.",
		call.=FALSE)
	distance
}

