smedian.sample <-
function(pos, v)
{
	w<-v[pos[1]:pos[2]][!is.na(v[pos[1]:pos[2]])]
	return(median(sample(w,length(w),replace=T),na.rm=T))
}
