get.chr.markers <-
function(chr.pair, marker.chr, paired.markers){

	#find where marker1 is on chromosome1 and marker2 is on chromosome2
	chr1.locale <- intersect(which(marker.chr[,1] == chr.pair[1]), which(marker.chr[,2] == chr.pair[2]))
	#and where marker1 is on chromosome2 and marker2 is on chromosome1
	chr2.locale <- intersect(which(marker.chr[,1] == chr.pair[2]), which(marker.chr[,2] == chr.pair[1]))
	
	#collect all these locations into a single vecto
	chr.pair.locale <- unique(c(chr1.locale, chr2.locale))
	chr.markers <- paired.markers[chr.pair.locale,,drop=FALSE]
	return(chr.markers)
	
	
}
