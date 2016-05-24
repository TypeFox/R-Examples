qsim <-
function (p, mu, sig){
	epi = 10^{-6}
	ll <- length(p)
	qq <- rep(0,ll)
	for (i in 1:ll){
		if (sig > 200){qq[i] <- qsim(p[i], mu, 200)}
		else{
		if (sig < 0.1){
			qq[i] <- qsim.norm(p[i], mu, sig)
			}
		#else if (sig > 100)
		else {	
		if ((sig<1 & mu*sig<0.1) | (sig>=1 & mu*sig <0.01)) {
			#norquan <- qsim.norm(p[i], mu, sig)
			grid <- c(seq(qsim.norm(p[i]/100, mu, sig), qsim.norm(min(p[i]+0.1, 1-epi), mu, sig), length.out=20), seq(0, mu, length.out=10)[2:10], seq(mu, 1, length.out=10)[2:9])
			} else if ((sig<1 & (1-mu)*sig<0.1) | (sig>=1 & (1-mu)*sig <0.01)) {
				#norquan <- qsim.norm(p[i], mu, sig)
				grid <- c(seq(qsim.norm(p[i]/100, mu, sig), qsim.norm(min(p[i]+0.1, 1-epi), mu, sig), length.out=20), seq(0, mu, length.out=10)[2:10], seq(mu, 1, length.out=10)[2:9])
				}
				else {
		step <- 0.1
		if (p[i] <= step/5) {grid <-  c(seq(10^(-8), epi, length.out=100), seq(0, mu, length.out=10)[2:10], seq(mu, 1, length.out=10)[2:9])}
		else if (p[i] >= (1-step/5)) {grid <-  c(seq(0, mu, length.out=10)[2:10], seq(mu, 1, length.out=10)[2:9], seq((1-epi), (1-10^(-8)),length.out=100))}
		else {grid <-  c(seq(0, mu, length.out=10)[2:10], seq(mu, 1, length.out=10)[2:9])}
		}
		ppgrid <- psim(grid, mu, sig)
		diffgrid <-sort(abs(ppgrid-p[i]), index.return=T)
		index <- diffgrid$ix
		sgrid <- grid[index]
		minidiff <- diffgrid$x[1]
		j <- 1
		while (minidiff>epi) {
			j <- j+1
			if (sgrid[1]==grid[1]){
				if (ppgrid[index[1]] > p[i]){
					lower <- sgrid[1]/2
					upper <- sgrid[1]
					}
			   else {
			   		lower <- min(sgrid[1:2])
			   		upper <- max(sgrid[1:2])
			   		}
				} else if (sgrid[1]==grid[length(grid)]) {
					if (ppgrid[index[1]] < p[i]) {
						lower <- sgrid[1] 
						upper <- (1+sgrid[1])/2
						}
					else {
						lower <- min(sgrid[1:2])
						upper <- max(sgrid[1:2])
						}
					} else {
						lower <- min(sgrid[1:2])
						upper <- max(sgrid[1:2])
						}
			step <- (upper-lower)/10
			grid <- seq(lower, upper, length.out=11)
			ppgrid <- psim(grid, mu, sig)
			diffgrid <- sort(abs(ppgrid-p[i]), index.return=T)
			index <- diffgrid$ix
			minidiff <- diffgrid$x[1]
			sgrid <- grid[index]
			if (j>100){break}
			}
		qq[i] <- sgrid[1]
		}}
		}
	return(qq)
	}
