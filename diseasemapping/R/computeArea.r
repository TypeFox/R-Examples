#return a list of areas
`computeArea`<-function(sp){
	lapply(sp@polygons, function(x) sum(sapply(x@Polygons, function(y) y@area * ifelse(y@hole, 0, 1))))
}
