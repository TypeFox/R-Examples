"mama" <-
function(dat, spl = TRUE) {
	dat <- data.frame(dat)
	plot <- as.character(dat[,1])
	spec <- as.character(dat[,2])
	if(ncol(dat) < 3){
		if(ncol(dat) < 2)
			stop("Mama needs at least two columns (vectors)!")
		dat$pres <- 1	
		}
	pres <- dat[,3]
	dat <- data.frame(plot, spec, pres)
	wide <- reshape(dat, v.names="pres", idvar="plot", timevar="spec", direction="wide")
	wide.nms <- sub("pres\\.", "", names(wide))
	if(spl){
		if(is.factor(pres)){
			wide <- sapply(c(1:ncol(wide)), function(x) as.character(wide[,x]))
		}
		wide[is.na(wide)] <- 0
	}
	rownames(wide) <- wide[,1]
	wide <- data.frame(wide)
	names(wide) <- wide.nms
	wide <- data.frame(wide[,-1])
	wide <- wide[order(rownames(wide)), ]
	wide <- wide[,order(names(wide))]
	wide
}