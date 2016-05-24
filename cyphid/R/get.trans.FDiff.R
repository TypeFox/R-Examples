get.trans.FDiff <-
function(dat, hz){
 
  set1 <- append(rev(dat[1:15]), dat)
  set2 <- append(set1, rev(set1[(length(set1)-14):length(set1)]))
  d1 <- rep(NA, length(set2))
  end2 <- length(set2) - 1

	for(h in 2:end2){
		d1[h] <- (set2[h+1]-set2[h-1]) / ((h+1) - (h-1))
		}
	
  
	af <- BWfilter(d1, hz,60) # 8
	
	d2 <- rep(NA, length(d1))
	
		for(k in 2:end2){
		d2[k] <- (af[k+1]-af[k-1]) / ((k+1) - (k-1))
	
	}
	amat2 <- d2
	CycleEnd <- length(amat2)-15
	acc <- amat2[16:CycleEnd]
	frames <- seq(1,length(dat))
  minindex <- get.peaks(-acc, 3)
	minframes <- frames[minindex]
	mindat <- acc[minindex]
	transmin <- which.min(mindat)
	return <- list(transmin=transmin, minframes=minframes, acc=acc)
	}

