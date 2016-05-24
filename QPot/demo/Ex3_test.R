# Below is a pipe every 10 characters.  The R Journal suggests that code is wrapped around 80 characters
#*********|*********|*********|*********|*********|*********|*********|*********

##### preparation #####
	# 0.0.1 write a function to create a legend for contour plots
	legend.col <- function(col, lev, xadj){ 
		opar <- par
		n <- length(col)
		bx <- par("usr")
		box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000 + 0, bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
		box.cy <- c(bx[3], bx[3])
		box.sy <- (bx[4] - bx[3]) / n
		xx <- rep(box.cx, each = 2) + xadj
		par(xpd = TRUE)
		for(i in 1:n){
			yy <- c(box.cy[1] + (box.sy * (i - 1)),
			box.cy[1] + (box.sy * (i)),
			box.cy[1] + (box.sy * (i)),
			box.cy[1] + (box.sy * (i - 1)))
			polygon(xx, yy, col = col[i], border = col[i])
		}
		par(new = TRUE)
		plot(0, 0, type = "n", ylim = c(min(lev), max(lev)), yaxt = "n", ylab = "",
		xaxt = "n", xlab = "", frame.plot = FALSE)
		axis(side = 4, las = 2, tick = FALSE, line = (.25 + xadj))
		par <- opar
	}
	# 0.0.2 preset color palettes for controus
	tsdens.col <- c("lightsteelblue", "white", "indianred")
	qp.col <- c("#FDE725FF", "#E3E418FF", "#C7E020FF", "#ABDC32FF", "#8FD744FF", "#75D054FF", "#5DC963FF", "#47C06FFF", "#35B779FF", "#28AE80FF", "#20A486FF", "#1F9A8AFF", "#21908CFF", "#24868EFF", "#287C8EFF", "#2C728EFF", "#31688EFF", "#355D8DFF", "#3B528BFF", "#404688FF", "#443A83FF", "#472D7BFF", "#481F71FF","#471163FF", "#440154FF")

##### 0.1 preperation #####
	# clean up and set seed
	rm(list=ls())
	your.favourite.number <- 3818919 #chris
	set.seed(your.favourite.number)

	#load libraries
	library(QPot)

##### 0.2 stochasticbounds.x simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "x*((1+alpha1)-x*x-x*y-y*y)"
		var.eqn.y <- "y*((1+alpha2)-x*x-x*y-y*y)"

	# 0.2.1 parameters
		model.state <- c(x = 0.5, y = 0.5)
		model.parms <- c(alpha1 = 1.25, alpha2 = 2)
		model.sigma <- 0.8
		model.time <- 5000
		model.deltat <- 0.01

	# 0.2.2 time series
	ts.ex3 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		temp.default.par <- par()
		TSPlot(ts.ex3, deltat = model.deltat)
		par(temp.default.par)
		TSPlot(ts.ex3, deltat = model.deltat, dim = 2 , line.alpha = 5)
		TSDensity(ts.ex3, dim = 1)
		TSDensity(ts.ex3, dim = 2 , contour.levels = 20 , contour.lwd = 0.1)

		# plots for the paper figures
			# Ex3_TS_1D.png
				# png(file = paste(plotwd, "/Ex3_TS_1D.png", sep = ""), width = 600, height = 350)
				# TSPlot(ts.ex3, deltat = model.deltat)
				# dev.off()
			# Ex3_TS_2D.png
				# png(file = paste(plotwd, "/Ex3_TS_2D.png", sep = ""), width = 400, height = 400)
				# par(mar = rep(2, 4), oma = rep(3,4))
				# TSPlot(ts.ex3, deltat = model.deltat, dim = 2, xlab ="", ylab = "", xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), line.alpha = 25)
				# mtext(expression(italic(x)), side = 1, line = 2.5)
				# mtext(expression(italic(y)), side = 2, line = 2.5)
				# dev.off()
			# Ex3_Dens_2D.png
				# k2 <- MASS::kde2d(ts.ex3[,2], ts.ex3[,3], n = 200)
				# k2dns <- k2$z/sum(k2$z)
				# k2cut <- cut(k2dns, 100, label = FALSE)
				# crramp <- colorRampPalette(tsdens.col)
				# colr <- crramp(100)
				# png(file = paste(plotwd, "/Ex3_Dens_2D.png", sep = ""), width = 400, height = 400)
				# par(mar = c(4, 4, 5, 5))
				# TSDensity(ts.ex3, dim = 2, xlab  ="", ylab = "", xlim = c(-3, 3), ylim = c(-3, 3), contour.levels = 20, contour.lines = T, las = 1, col2d = tsdens.col , contour.lwd = 0.1, kde2d.n = 200, xaxs = "i", yaxs = "i")
				# legend.col(col = colr, lev = k2dns, xadj = 0.1)
    				# mtext(expression(italic(x)), side = 1, line = 2.5)
				# mtext(expression(italic(y)), side = 2, line = 2.5)
				# dev.off()


##### 0.3 local quasi-potential!!! #####
		equation.x = "x*((1+1.25)-x*x-x*y-y*y)"
		equation.y = "y*((1+2)-x*x-x*y-y*y)"
		bounds.x = c(-3, 3)
		bounds.y = c(-3, 3)
		step.number.x = 6000
		step.number.y = 6000
		eq1.x = 0
		eq1.y = -1.73205
		eq2.x = 0
		eq2.y = 1.73205

eq1.local <- QPotential(x.rhs = equation.x, x.start = eq1.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq1.y, y.bound = bounds.y, y.num.steps = step.number.y)
eq2.local <- QPotential(x.rhs = equation.x, x.start = eq2.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq2.y, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
ex3.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local),unstable.eq.x = c(0, -1.5, 1.5), unstable.eq.y = c(0, 0, 0), x.bound = bounds.x, y.bound = bounds.y)


##### 0.5 quasi-potential vizualization!!! #####
QPContour(ex3.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 1)

	# figure for the paper
		# png(file = paste(plotwd, "/Ex3_QP_contour.png", sep = ""), width = 400, height = 400)
		# par(oma = c(0,1,0,2))
		# QPContour(surface = ex3.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 4, xlab = expression(italic("x")), ylab = expression(italic("y")), n.contour.lines = 20)
		# k2dns <- seq(min(ex3.global, na.rm = T), max(ex3.global, na.rm = T), 0.01)
		# k2cut <- cut(k2dns, 100, label = FALSE)
		# crramp <- colorRampPalette(qp.col)
		# colr <- crramp(100)
		# legend.col(col = colr, lev = k2dns, xadj = 0.1)
		# dev.off()

##### 0.6 vector field decompisition!!! #####
	# 0.6.0 all fields
	VDAll <- VecDecomAll(surface = ex3.global, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
	VecDecomPlot(x.field = VDAll[,,1], y.field = VDAll[,,2], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", tail.length = 0.2, head.length = 0.03)
	VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 1/3, head.length = 0.03)
	VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.5, head.length = 0.03)

	# 0.6.1 vector field
	VDV <- VecDecomVec(x.num.steps = step.number.y, y.num.steps= step.number.y, x.rhs=equation.x, y.rhs=equation.y, x.bound= bounds.x, y.bound= bounds.x)
		VecDecomPlot(x.field=VDV[,,1], y.field=VDV[,,2], dens=c(50,50), x.bound= bounds.x, y.bound=bounds.y , tail.length = 0.2 , length = 0.03)

	# 0.6.2 gradient field	
	VDG <- VecDecomGrad(e3.global)
		VecDecomPlot(x.field=VDG[,,1], y.field=VDG[,,2], dens=c(25,25), x.bound= bounds.x, y.bound= bounds.y , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")

	# 0.6.3 remainder field
	VDR <- VecDecomRem(surface=e3.global, x.rhs=equation.x, y.rhs=equation.y, x.bound=bounds.x, y.bound=bounds.y)
		VecDecomPlot(x.field=VDR[,,1], y.field=VDR[,,2], dens=c(25,25), x.bound= bounds.y, y.bound=bounds.y , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")

##### 0.7 3D graphs!!! #####
	library(rgl)
	dens.sub <- c(1000,1000)
	global.sub <- ex3.global[round(seq(1,nrow(ex3.global),length.out=dens.sub[1])) , round(seq(1,ncol(ex3.global),length.out=dens.sub[2]))]
	# global.sub[global.sub > 2] <- NA # to cut off large values for a better view of the basin
	persp3d(x = global.sub, col = "orange", expand = 1.1)
