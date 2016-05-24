## script for transforming from Braun-Blanquet to
## numeric data...
"bb2num" <-
function(dat, from = c("r", "+", "1", "2", "3", "4", "5"), to = c(0.1, 1, 5, 15, 37.5, 62.5, 87.5)){
	odd <- function(x){x%%2 == 1}
	if(sum((length(from)!=length(to)), (length(from)!=2*length(to))) == 2){
		stop("from and to have to have the same length")
		}
	r.nms <- rownames(dat)
	c.nms <- names(dat)
	class.from <- eval(parse(text=paste("as", class(from), sep=".")))
	dat.new <- apply(dat,  2, class.from)
	if(is.character(from)){
		for(i in 1:length(from)){
			dat.new[dat==from[i]] <- to[i]
		}
	}
	else{
		vek <- c(1:length(from))
		to <- as.vector(sapply(to, rep, 2))
		for(i in vek[odd(vek)]){
			dat.new[(dat > from[i]) & (dat <= from[i+1])] <- to[i]
		}
	}
	class.to <- eval(parse(text=paste("as", class(to), sep=".")))
	dat <- apply(dat.new, 2, class.to)
	dat[is.na(dat)] <- 0
	dat <- data.frame(dat)
	dimnames(dat) <- list(r.nms, c.nms)
	dat
	}