smedmad <-
function(pos,v)
	c(median(v[pos[1]:pos[2]],na.rm=T),mad(v[pos[1]:pos[2]],na.rm=T))
