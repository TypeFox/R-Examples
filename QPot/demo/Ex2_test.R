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


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "-(y-beta) + mu*(x-alpha)*(1-(x-alpha)^2-(y-beta)^2) "
		var.eqn.y <- "(x-alpha) + mu*(y-beta)*(1-(x-alpha)^2-(y-beta)^2)"

	# 0.2.1 parameters
		model.state <- c(x = 3, y = 3)
		model.parms <- c(alpha = 4, beta = 5, mu = 0.2)
		model.sigma <- 0.1
		model.time <- 2500
		model.deltat <- 0.005

	# 0.2.2 time series
ts.ex2 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		temp.default.par <- par()
		TSPlot(ts.ex2, deltat = model.deltat)
		par(temp.default.par)
		TSPlot(ts.ex2, deltat = model.deltat, dim = 2, line.alpha = 25)
		TSDensity(ts.ex2, dim = 1)
		TSDensity(ts.ex2, dim = 2)

		# plots for the paper figures
			print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
			# Ex2_TS_1D.png
				# png(file = paste(plotwd, "/Ex2_TS_1D.png", sep = ""), width = 600, height = 350)
				# TSPlot(ts.ex2, deltat = model.deltat)
				# dev.off()
			# Ex2_TS_2D.png
				# png(file = paste(plotwd, "/Ex2_TS_2D.png", sep = ""), width = 400, height = 400)
				# par(mar = rep(2, 4), oma = rep(3,4))
				# TSPlot(ts.ex2, deltat = model.deltat, dim = 2, xlab ="", ylab = "", xlim = c(2.5, 6.5), ylim = c(2.5, 6.5), line.alpha = 25)
				# mtext(expression(italic(x)), side = 1, line = 2.5)
				# mtext(expression(italic(y)), side = 2, line = 2.5)
				# dev.off()
			# Ex2_Dens_2D.png
				# k2 <- MASS::kde2d(ts.ex2[,2], ts.ex2[,3], n = 200)
				# k2dns <- k2$z/sum(k2$z)
				# k2cut <- cut(k2dns, 100, label = FALSE)
				# crramp <- colorRampPalette(tsdens.col)
				# colr <- crramp(100)

				# png(file = paste(plotwd, "/Ex2_Dens_2D.png", sep = ""), width = 400, height = 400)
				# par(mar = c(4, 4, 5 , 5))
				# TSDensity(ts.ex2, dim = 2, xlab  ="", ylab = "", xlim = c(2.5, 6.5), ylim = c(2.5, 6.5), contour.levels = 25, contour.lines = T, las = 1, col2d = tsdens.col, contour.lwd = 0.25, kde2d.n = 200, xaxs = "i", yaxs = "i")
				# legend.col(col = colr, lev = k2dns, xadj = 0.1)
    				# mtext(expression(italic(x)), side = 1, line = 2.5)
				# mtext(expression(italic(y)), side = 2, line = 2.5)
				# dev.off()

##### 0.3 local quasi-potential!!! #####
		equation.x = "-(y-5) + 0.2*(x-4)*(1-(x-4)^2-(y-5)^2)"
		equation.y = "(x-4) + 0.2*(y-5)*(1-(x-4)^2-(y-5)^2)"
		bounds.x = c(-0.5, 7.5)
		bounds.y = c(-0.5, 7.5)
		step.number.x = 4000
		step.number.y = 4000
		xinit = 4.15611
		yinit = 5.987774

eq1.qp <- QPotential(x.rhs = equation.x, x.start = xinit, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = yinit, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
	#same as the local quasi-potential, calculated above


##### 0.5 quasi-potential vizualization!!! #####
QPContour(eq1.qp, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 10, xlab=expression(italic(x)), ylab=expression(italic(y)))
	# figure for the paper
		# png(file = paste(plotwd, "/Ex2_QP_contour.png", sep = ""), width = 400, height = 400)
		# par(oma = c(0,1,0,2))
		# QPContour(surface = eq1.qp, dens = c(1000, 1000), x.bound = c(-0.5, 7.5), y.bound = c(-0.5, 7.5), c.parm = 10, xlab = expression(italic("x")), ylab = expression(italic("y")), col.contour = qp.col, n.contour.lines = 15, xlim = c(1.5, 7.5), ylim = c(1.5, 7.5))
		# k2dns <- seq(min(eq1.qp, na.rm = T), max(eq1.qp, na.rm = T), 0.01)
		# k2cut <- cut(k2dns, 100, label = FALSE)
		# crramp <- colorRampPalette(qp.col)
		# colr <- crramp(100)
		# legend.col(col = colr, lev = k2dns, xadj = 0.1)
		# dev.off()


##### 0.6 vector field decompisition!!! #####
	# 0.6.1 vector field
	VDV <- VecDecomVec(x.num.steps = step.number.x, y.num.steps = step.number.y, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
		VecDecomPlot(x.field = VDV[,,1], y.field = VDV[,,2], dens = c(25,25), x.bound = bounds.x, y.bound = bounds.y, tail.length = 0.5, head.length = 0.03, arrow.type = "proportional", xlim = c(2, 6), y.lim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))

	# 0.6.2 gradient field	
	VDG <- VecDecomGrad(eq1.qp)
		VecDecomPlot(x.field = VDG[,,1], y.field = VDG[,,2], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, tail.length = 0.2, head.length = 0.05 , arrow.type = "proportional", xlim = c(2, 6), ylim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))

	# 0.6.3 remainder field
	VDR <- VecDecomRem(surface = eq1.qp, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.x)
		VecDecomPlot(x.field = VDR[,,1], y.field = VDR[,,2], dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.x , tail.length = 0.2 , head.length = 0.03 , arrow.type = "proportional", xlim = c(2, 6), ylim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))
