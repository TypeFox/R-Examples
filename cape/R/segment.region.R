segment.region <-
function(region.min, region.max, num.points, alignment = c("center", "ends")){
	
	if(num.points < 2){
		return(mean(c(region.min, region.max)))
		}
	
	if(length(grep("c", alignment)) > 0){
		alignment <- "center"
		}
	
	
	total.region <- region.max - region.min
	
	if(alignment == "ends"){
		point.seq <- seq(region.min, region.max, total.region/(num.points-1))
		return(point.seq)
		}


	if(alignment == "center"){
		#first break the segment into n+1 regions
		point.seq <- seq(region.min, region.max, total.region/num.points)
		#find the center of each region
		cons.pairs <- consec.pairs(1:length(point.seq))
		center.points <- apply(cons.pairs, 1, function(x) mean(c(point.seq[x[1]], point.seq[x[2]])))
		return(center.points)
		}
		
		
	# plot(center.points, rep(1, length(center.points)), xlim = c(region.min, region.max), col = "red")
	# points(point.seq, rep(1.2, length(point.seq)), col = "blue")
	
	
}
