#
# A derivate function 
# 
deriva <- function(b,funcpa){
	m <- length(b)
	bh2 <- bh <- rep(0,m)
	v <- rep(0,(m*(m+3)/2))
	fcith <- fcith2 <- rep(0,m)
	# function 
	rl <- funcpa(b)
	# gradient null
	for(i in 1:m){
		bh <- bh2 <- b
		th <- max(1E-7,(1E-4*abs(b[i])))
		bh[i] <- bh[i] + th
		bh2[i] <- bh2[i] - th
	
		fcith[i] <- funcpa(bh)
		fcith2[i] <- funcpa(bh2)	
	}
	k <- 0
	m1 <- m*(m+1)/2
	l <- m1
	for(i in 1:m){
		l <- l+1
		bh <- b
		thn <- - max(1E-7,(1E-4*abs(b[i])))
		v[l] <- -(fcith[i]-fcith2[i])/(2*thn)  
	
		for(j in 1:i){
			bh <- b
			k <- k+1
			thi <- max(1E-7,(1E-4*abs(b[i])))
			thj <- max(1E-7,(1E-4*abs(b[j])))
			th <- thi * thj
			bh[i] <- bh[i]+thi
			bh[j] <- bh[j]+thj
			temp <-funcpa(bh)
			v[k] <- -(temp-(fcith[j])-(fcith[i])+rl)/th
			
		}
	}	
	result <- list(v=v,rl=rl)
	return(result)
}

