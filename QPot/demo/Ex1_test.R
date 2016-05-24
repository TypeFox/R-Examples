# Below is a pipe every 10 characters.  The R Journal suggests that code is wrapped around 80 characters
#*********|*********|*********|*********|*********|*********|*********|*********
require(QPot)
require(rgl)

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
	# 0.1.0 clean up and set seed
	rm(list=ls())
	your.favourite.number <- 3818919 #chris
	set.seed(your.favourite.number)


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa+(x^2)))"
		var.eqn.y <- "((gamma*(x^2)*y)/(kappa+(x^2))) - mu*(y^2)"

	# 0.2.1 parameters
		model.state <- c(x = 1, y = 2)
		model.parms <- c(alpha = 1.54, beta = 10.14, delta = 1, gamma = 0.476, kappa = 1, mu = 0.112509)
		model.sigma <- 0.05
		model.time <- 12500
		model.deltat <- 0.025

	# 0.2.2 time series
	ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		temp.default.par <- par()
		TSPlot(ts.ex1, deltat = model.deltat)
		par(temp.default.par)
		TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
		TSDensity(ts.ex1, dim = 1)
		TSDensity(ts.ex1, dim = 2)

		# plots for the paper figures
		# time series, 1D
			# png(file = paste(plotwd, "/Ex1_TS_1D.png", sep = ""), width = 600, height = 350)
			# TSPlot(ts.ex1, deltat = model.deltat)
			# dev.off()
		# time series, 2D
			# png(file = paste(plotwd, "/Ex1_TS_2D.png", sep = ""), width = 400, height = 400)
			# par(mar = rep(2, 4), oma = rep(3,4))
			# TSPlot(ts.ex1, deltat = model.deltat, dim = 2, xlab ="", ylab = "", xlim = c(0.5, 6), ylim = c(0.5, 6), line.alpha = 50)
			# mtext(expression(italic(x)), side = 1, line = 2.5)
			# mtext(expression(italic(y)), side = 2, line = 2.5)
			# dev.off()
		# density, 2D
			# k2 <- MASS::kde2d(ts.ex1[,2], ts.ex1[,3], n = 200)
			# k2dns <- k2$z/sum(k2$z)
			# k2cut <- cut(k2dns, 100, label = FALSE)
			# crramp <- colorRampPalette(tsdens.col)
			# colr <- crramp(100)
			# png(file = paste(plotwd, "/Ex1_Dens_2D.png", sep = ""), width = 400, height = 400)
			# par(mar = c(4, 4, 4 , 4))
			# TSDensity(ts.ex1, dim = 2, xlab  ="", ylab = "", xlim = c(0.5, 6), ylim = c(0.5, 6), contour.levels = 25, contour.lines = T, las = 1, col2d = tsdens.col, contour.lwd = 0.2, kde2d.n = 200, xaxs = "i", yaxs = "i")
			# legend.col(col = colr, lev = k2dns, xadj = 0.1)
    			# mtext(expression(italic(x)), side = 1, line = 2.5)
			# mtext(expression(italic(y)), side = 2, line = 2.5)
			# dev.off()


##### 0.3 local quasi-potential!!! #####
	equation.x <- Model2String(var.eqn.x, parms = model.parms)
	# Could also input the values by hand and use this version
	# equation.x = "1.54*x*(1.0-(x/10.14))-(y*(x^2))/(1.0+(x^2))"
	equation.y <- Model2String(var.eqn.y, parms = model.parms,
	supress.print =TRUE) # does not print to screen
	# equation.y = "((0.476*(x^2)*y)/(1+(x^2)))-0.112509*(y^2)"
	bounds.x = c(-0.5, 20.0)
	bounds.y = c(-0.5, 20.0)
	step.number.x = 4100
	step.number.y = 4100
	eq1.x = 1.40491
	eq1.y = 2.80808
	eq2.x = 4.9040
	eq2.y = 4.06187

	eq1.local <- QPotential(x.rhs = equation.x, x.start = eq1.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq1.y,  y.bound = bounds.y, y.num.steps = step.number.y)
	eq2.local <- QPotential(x.rhs = equation.x, x.start = eq2.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq2.y, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
	ex1.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local), unstable.eq.x = c(0, 4.2008), unstable.eq.y = c(0, 4.0039), x.bound = bounds.x, y.bound = bounds.y)


##### 0.5 quasi-potential vizualization!!! #####
	QPContour(surface = ex1.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5)
	# A paper-like graph
		# par(mfrow = c(1, 2), mar = c(4, 2, 0.5, 0.5), oma = rep(3,4))
		# k2dns <- seq(min(ex1.global, na.rm = T), max(ex1.global, na.rm = T), 0.001)
		# k2cut <- cut(k2dns, 100, label = FALSE)
		# crramp <- colorRampPalette(qp.col)
		# colr <- crramp(100)

		# png(file = paste(plotwd, "/Ex1_QP_contour.png", sep = ""), width = 800, height = 400)
		# par(mfrow = c(1, 2), mar = c(4, 0, 2, 2), oma = c(0, 4, 2, 2))
		# QPContour(surface = ex1.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 1, xlab = expression(italic("x")), ylab = "")
		# mtext(text = expression(italic("y")), side = 2, line = 2.5, outer = T, at = 0.5)
		# QPContour(surface = ex1.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5, xlab = expression(italic("x")), ylab = "")
		# legend.col(col = colr, lev = k2dns, xadj = 0.1)
		# dev.off()


##### 0.6 vector field decompisition!!! #####
	# 0.6.0 all fields
	VDAll <- VecDecomAll(surface = ex1.global, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
	VecDecomPlot(x.field = VDAll[,,1], y.field = VDAll[,,2], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 11), ylim = c(0, 6), arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
	VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.25, head.length = 0.025)
	VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.35, head.length = 0.025)

	# 0.6.1 vector field
	VDV <- VecDecomVec(x.num.steps = step.number.y, y.num.steps= step.number.y, x.rhs=equation.x, y.rhs=equation.y, x.bound= bounds.x, y.bound= bounds.x)
	VecDecomPlot(x.field=VDV[,,1], y.field=VDV[,,2], dens=c(25,25), x.bound= bounds.x, y.bound=bounds.y,
		xlim = c(0, 11), ylim = c(0, 6), arrow.type = "equal", tail.length = 0.25, head.length = 0.025)

	# 0.6.2 gradient field	
	VDG <- VecDecomGrad(ex1.global)
		VecDecomPlot(x.field=VDG[,,1], y.field=VDG[,,2], dens=c(25,25), x.bound= bounds.x, y.bound= bounds.y , tail.length = 0.25, head.length = 0.025, arrow.type = "proportional")

	# 0.6.3 remainder field
	VDR <- VecDecomRem(surface=ex1.global, x.rhs=equation.x, y.rhs=equation.y, x.bound=bounds.x, y.bound=bounds.y)
		VecDecomPlot(x.field=VDR[,,1], y.field=VDR[,,2], dens=c(25,25), x.bound= bounds.y, y.bound=bounds.y , tail.length = 0.35, head.length = 0.025, arrow.type = "proportional")

	# plots for the paper figure
		# Paper vector decomposition plot
			# png(file = paste(plotwd, "/Ex1_VecDecom.png", sep = ""), width = 500, height = 500)
			# par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(3,3,1,1))
			# VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.25, head.length = 0.025)
			# mtext(side = 3, text = expression(Gradient~field~(-nabla~phi(x, y))))
			# mtext(side = 2, text = "Proportional arrow lengths", line = 3.75)
			# mtext(side = 2, text = expression(italic(y)), line = 2.25)
			# VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.35, head.length = 0.025)
			# mtext(side = 3, text = expression(Remainder~field~(bold(r)(x, y))))
			# VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", tail.length = 0.15, head.length = 0.025)
			# mtext(side = 2, text = "Equal arrow lengths", line = 3.75)
			# mtext(side = 2, text = expression(italic(y)), line = 2.25)
			# mtext(side = 1, text = expression(italic(x)), line = 2.25)	
			# VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", tail.length = 0.15, head.length = 0.025)
			# mtext(side = 1, text = expression(italic(x)), line = 2.25)
			# dev.off()



##### 0.7 3D graphs!!! #####
	library(rgl)
	dens.sub <- c(4000,4000)
	global.sub <- ex1.global[round(seq(1,nrow(ex1.global),length.out=dens.sub[1])) , round(seq(1,ncol(ex1.global),length.out=dens.sub[2]))]
	# global.sub[global.sub > 0.02] <- NA # to cut off large values for a better view of the basin
	persp3d(x = global.sub, col = "orange", expand = 1.1, xlim = c(0.05, 0.35), ylim = c(0.1, 0.3), zlim = c(0, 0.01), xlab = "X", ylab = "Y", zlab = intToUtf8(0x03A6))
	# rgl.snapshot(filename = paste(plotwd, "/Ex1_3D.png", sep = ""), fmt = "png", top = TRUE ) #to print view as a png
