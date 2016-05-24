
calc.props <- function(y){

	sums <- apply(y, 1, sum)
	sums <- ifelse(sums == 0, 1, sums)
	p <- sweep(y, 1, FUN = "/", sums)

	return(p)

}


calc.color <- function(colors, levels, x){

	z <- ifelse(x == 0, 1, which(x <= levels)[1]-1)

	return(colors[z])

}



click.plot <- function(X, y = NULL, file = NULL, id, states = NULL, marg = 1, font.cex = 2, font.col = "black", cell.cex = 1, cell.lwd = 1.3, cell.col = "black", sep.lwd = 1.3, sep.col = "black", obs.lwd = NULL, colors = c("lightcyan", "pink", "darkred"), col.levels = 8, legend = TRUE, leg.cex = 1.3, top.srt = 0, frame = TRUE){

	K <- max(id)
	p <- dim(X)[1]
	n <- dim(X)[3]

	P1 <- X

	for (i in 1:n) P1[,,i] <- calc.props(X[,,i])

	Nk <- NULL

	last <- 0
	if (!is.null(y)) y.new <- NULL


	if (is.null(obs.lwd)){ # calculate median colors

		P.med <- array(rep(NA, p^2*K), c(p, p, K))		

		for (k in 1:K){

			ind <- which(id == k)
			nk <- length(ind)
			Nk <- c(Nk, nk)
			P.med[,,k] <- apply(P1[,,ind], c(1,2), median)

		}

	} else { 	# sort data according to group assignments

		P <- P1		

		for (k in 1:K){

			ind <- which(id == k)
			nk <- length(ind)
			Nk <- c(Nk, nk)
			P[,,(last+1):(last+nk)] <- P1[,,ind]
			last <- last + nk

			if (!is.null(y)){
				y.new <- c(y.new, y[ind])
			}

		}

	}

	Nk.cum <- cumsum(Nk)


	colors <- colorRampPalette(colors)(col.levels)

	levels = seq(0.0, 1.0, length.out = col.levels+1)

	grid <- seq(0, 1, length.out = p + 1)
	grid.step <- grid[2] - grid[1]

	par(mar = rep(0.1, 4))
	if (legend){
		if (is.null(y)){
			plot( c(-grid.step/2 * marg, 1), c(-grid.step, 1 + grid.step/2 * marg), type = "n", xlab = "", ylab = "", axes = FALSE)
		} else {
			plot( c(-grid.step/2 * marg, 1 + grid.step), c(-grid.step, 1 + grid.step/2 * marg), type = "n", xlab = "", ylab = "", axes = FALSE)
		}
	} else {
		if (is.null(y)){
			plot( c(-grid.step/2 * marg, 1), c(0, 1 + grid.step/2 * marg), type = "n", xlab = "", ylab = "", axes = FALSE)
		} else {
			plot( c(-grid.step/2 * marg, 1 + grid.step), c(0, 1 + grid.step/2  * marg), type = "n", xlab = "", ylab = "", axes = FALSE)
		}
	}
	if (frame) box()

	# state numbers

	y1 <- 1 + grid.step / 3 * marg
	for (j in 1:p){
		x1 <- (grid[j] + grid[j+1]) / 2
		if (is.null(states)){
			text(x1, y1, j, cex = font.cex, col = font.col, srt = top.srt)
		} else {
			text(x1, y1, states[j], cex = font.cex, col = font.col, srt = top.srt)
		}
	}

	x1 <- -grid.step / 3 * marg
	for (j in 1:p){
		y1 <- (grid[j] + grid[j+1]) / 2
		if (is.null(states)){
			text(x1, y1, p+1-j, cex = font.cex, col = font.col)
		} else {
			text(x1, y1, states[p+1-j], cex = font.cex, col = font.col)
		}
	}

	# margin between cells
	eps <- grid.step / 20 / cell.cex






	if (is.null(obs.lwd)){	# median color polygons

		Nk.cum <- c(0, Nk.cum)
		step <- (grid.step - 2 * eps) / n

		for (j in 1:p){
			x1 <- grid[j]
			y1 <- grid[j]
			for (i in 1:p){
				for (k in 1:K){

					polygon(c(grid[j]+eps, grid[j+1]-eps, grid[j+1]-eps, grid[j]+eps),
						c(grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step,
						  grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step,
						  grid[p-i+1]+eps+(Nk.cum[k+1]+0.5)*step,
						  grid[p-i+1]+eps+(Nk.cum[k+1]+0.5)*step),
						  col = calc.color(colors, levels, P.med[i,j,k]),
						  border = sep.col, lwd = sep.lwd)

				}
			}
		}


		# cell frames

		for (j in 1:p){
			x1 <- grid[j]
			y1 <- grid[j]
			for (i in 1:p){
				polygon(c(grid[j]+eps, grid[j+1]-eps, grid[j+1]-eps, grid[j]+eps),
					c(grid[p-i+1]+eps, grid[p-i+1]+eps, grid[p-i+2]-eps, grid[p-i+2]-eps),
					border = cell.col, lwd = cell.lwd)
			}

		}


		if (!is.null(y)){	# additional column of cells to represent betas

			for (j in 1:p){
				for (k in 1:K){
					polygon(c(max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2,
						  max(grid) + grid.step / 2 + 5*eps*1.2, max(grid) + 5*eps*1.2),
						c(grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step,
						  grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step,
						  grid[p-j+1]+eps+(Nk.cum[k+1]+0.5)*step,
						  grid[p-j+1]+eps+(Nk.cum[k+1]+0.5)*step),
						col = calc.color(colors, levels, mean(y[id == k] == j)),
						border = sep.col, lwd = sep.lwd)
				}
			}


			# cell frames

			for (j in 1:p){
				polygon(c(max(grid) + grid.step / 2 + 5*eps*1.2, max(grid) + 5*eps*1.2,
					  max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2),
					c(grid[p-j+1]+eps, grid[p-j+1]+eps, grid[p-j+2]-eps, grid[p-j+2]-eps),
					  border = cell.col, lwd = cell.lwd)
			}

		}





	} else { 	# observation lines

		step <- (grid.step - 2 * eps) / n
		for (j in 1:p){
			for (i in 1:p){
				curr <- P[i,j,]
				for (h in 1:n){
					lines(c(grid[j]+eps*1.2, grid[j+1]-eps*1.2),
					      c(grid[p-i+1]+eps+h*step, grid[p-i+1]+eps+h*step),
					      col = calc.color(colors, levels, curr[h]), lwd = obs.lwd)
				}
				if (K != 1){
					for (k in 1:(K-1)){
						lines(c(grid[j]+eps*1.2, grid[j+1]-eps*1.2),
						      c(grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step,
							grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step),
						      col = sep.col, lwd = sep.lwd)
					}
				}
			}
		}


		# cell frames

		for (j in 1:p){
			x1 <- grid[j]
			y1 <- grid[j]
			for (i in 1:p){
				polygon(c(grid[j]+eps, grid[j+1]-eps, grid[j+1]-eps, grid[j]+eps),
					c(grid[p-i+1]+eps, grid[p-i+1]+eps, grid[p-i+2]-eps, grid[p-i+2]-eps),
					border = cell.col, lwd = cell.lwd)
			}

		}



		if (!is.null(y)){	# additional column of cells to represent betas

			for (j in 1:p){
				for (i in 1:n){
					lines(c(max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2),
					      c(grid[p-j+1]+eps+i*step, grid[p-j+1]+eps+i*step),
					      col = calc.color(colors, levels, y.new[i] == j), lwd = obs.lwd)
				}

				if (K != 1){
					for (k in 1:(K-1)){
						lines(c(max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2),
						      c(grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step,
							grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step),
					     		col = sep.col, lwd = sep.lwd)
					}	
				}
			}


			# cell frames

			for (j in 1:p){
				polygon(c(max(grid) + grid.step / 2 + 5*eps*1.2, max(grid) + 5*eps*1.2,
					  max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2),
					c(grid[p-j+1]+eps, grid[p-j+1]+eps, grid[p-j+2]-eps, grid[p-j+2]-eps),
					  border = cell.col, lwd = cell.lwd)
			}

		}



	}


	# construct legend

	if (legend){

		sep.X <- seq(0+eps, 1-eps, length.out = col.levels + 1)

		for (j in 1:col.levels){
			polygon(c(sep.X[j], sep.X[j+1], sep.X[j+1], sep.X[j]),
				c(-grid.step*0.5, -grid.step*0.5, -grid.step*0.25, -grid.step*0.25),
				border = cell.col, lwd = cell.lwd, col = colors[j])
			text(sep.X[j], -grid.step*0.75, round(levels[j], digits = 2), cex = leg.cex)
		}
		text(1-eps, -grid.step*0.75, 1, cex = leg.cex, col = font.col)

	}


	if (!is.null(file)) dev.copy2pdf(file = file)

}


