plot.fpt.density <-
function (x, from.t0, to.T, dp.legend = TRUE, dp.legend.cex = 1, ylab = TRUE, growth.points = FALSE, instants = FALSE, ...) 
{
    	if (!is.fpt.density(x)) 
        	stop(paste(sQuote("x"), "is not of class", shQuote("fpt.density")))
	
    	ncp <- !is.null(x$y.x0)
	
	Args <- as.list(attr(x, "Call"))
    	if (is.element("from.t0", names(Args))) fromt0Call <- Args$from.t0 else fromt0Call <- FALSE
    	if (missing(from.t0)) from.t0 <- fromt0Call
    	if (is.element("to.T", names(Args))) toTCall <- Args$to.T else toTCall <- FALSE
    	if (missing(to.T)) to.T <- toTCall

	if (ncp){
		Call <- attr(x, "Call")
		Args <- as.list(Call[c("t0","T", "id")])
	  	S <- Call$S
		if (is.list(Args$id)) Args$id <- Args$id[[2]] else Args$id <- attr(attr(x, "summary.fptl"), "id")[[2]]
	  	env <- Call$env
	  	if (!is.null(env)){
			if (is.call(env)) Args <- c(Args, eval(env))
			else Args <- c(Args, attr(attr(x, "summary.fptl"), "vars")[[as.character(env)]])	
	  	}		
		dp <- attr(attr(x, "summary.fptl"), "dp")
		Args$id <- eval(parse(text = paste("substitute(", Args$id, ")", sep = "")))
		t1 <- min(sapply(attr(x, "summary.fptl"), function(x) x$instants[1, 2]))
		t2 <- max(sapply(attr(x, "summary.fptl"), function(x) x$instants[length(x)]))

		if (growth.points) warning(paste(sQuote("growth.points"), "argument has been ignored."), call. = FALSE)
		if (instants) warning(paste(sQuote("instants"), "argument has been ignored."), call. = FALSE)
		growth.points <- FALSE
		instants <- FALSE
	}
	else{
		Args <- formals(FPTL)[c("t0","T", "x0")]
		Call <- attr(attr(x, "summary.fptl"), "FPTLCall")[[1]]
		Args <- as.list(Call[c("t0","T", "x0")])
	  	S <- Call$S
		env <- Call$env
	  	if (!is.null(env)){
			if (is.call(env)) Args <- c(Args, eval(env))
			else Args <- c(Args, attr(attr(x, "summary.fptl"), "vars")[[as.character(env)]])	
	  	}		
		dp <- attr(attr(x, "summary.fptl"), "dp")
		t1 <- attr(x, "summary.fptl")[[1]]$instants[1, 2]
		t2 <- attr(x, "summary.fptl")[[1]]$instants[length(attr(x, "summary.fptl")[[1]]$instants)]		
	}
	
	par(cex = 1, ps = 9)
    	par(mar = c(3 + 2 * instants, 3.25 + ylab, 2.5, 1) + 0.1, ...)        
    	pin.height <- par("pin")[2] - sum(par("mai")[c(1,3)]*(par("cex")-1))
    	if (dp.legend) {
		dp.labels <- vector("list", 3)
		if (ncp) dp.labels[[1]] <- parse(text = deparse(substitute(paste(" Diffusion process: ", list(group("{", list(X(t), paste(t, "  in  ", group("[", 
                							list(t0, T), "]"))), "}"), ~~X(t0) %~% id)), Args)))
		else dp.labels[[1]] <- parse(text = deparse(substitute(paste(" Diffusion process: ", list(group("{", list(X(t), paste(t, "  in  ", group("[", 
                					list(t0, T), "]"))), "}"), ~~P(X(t0) == x0) == 1)), Args)))
	  	logic <- sapply(Args, is.character)
	  	Args[logic] <- sapply(Args[logic], function(x) parse(text = x)[[1]])
		dpMean <- eval(parse(text = paste("substitute(", dp$mean, ", Args[logic])", sep = "")))
	  	dpVar <- eval(parse(text = paste("substitute(", dp$var, ", Args[logic])", sep = "")))
	  	S <- eval(parse(text = paste("substitute(", S, ", Args[logic])", sep = "")))
	  	logic <- sapply(Args, is.numeric) & (sapply(Args, length) == 1L)
	  	dpMean <- eval(parse(text = paste("substitute(", deparse(dpMean, width.cutoff=500L), ", Args[logic])", sep = "")))
	  	dpVar <- eval(parse(text = paste("substitute(", deparse(dpVar, width.cutoff=500L), ", Args[logic])", sep = "")))
	  	S <- eval(parse(text = paste("substitute(", deparse(S, width.cutoff=500L), ", Args[logic])", sep = "")))
	  	dp.labels[[2]] <- parse(text = gsub("*", "%.%", deparse(substitute(list(paste(A[1](x, t) == m, " "), ~~A[2](x, t) == v), list(m = dpMean, v = dpVar))), fixed=TRUE))					
        	dp.labels[[3]] <- parse(text = gsub("*", "%.%", deparse(substitute(paste(" Boundary: ", S(t) == s), list(s = S))), fixed=TRUE))		
	  	dp.width <- sum(sapply(dp.labels[1:2], strwidth, units = "inch", cex = dp.legend.cex)) + strwidth(" ,  ", units = "inch", cex = dp.legend.cex)
	  	pin.width <- par("pin")[1] - sum(par("mai")[c(2,4)]*(par("cex")-1))	  
	  	dp.h <- sapply(dp.labels, strheight, units = "inch", cex = dp.legend.cex)
	  	logic <- (dp.width < pin.width)
	  	if (logic) dp.height <- max(dp.h[1:2]) + dp.h[3] + 0.3*par("cin")[2] else dp.height <- sum(dp.h) + 0.4*par("cin")[2]	      
    	}
    	else dp.height <- numeric(1)

	gp.h <- 0.03 * pin.height
	ti.h <- gp.h
	if (growth.points) gp.h <- max(gp.h, 1.25 * strheight(expression(t), units = "inch", cex = sqrt(par("cex"))) + 0.2*par("cin")[2])
	if (instants) ti.h <- max(ti.h, 1.5 * strheight(expression(t[1]), units = "inch", cex = sqrt(par("cex"))) + 0.2*par("cin")[2])
    
	if ((!from.t0) & fromt0Call & (x$x[1] < t1)) {
		lg <- (x$x >= t1)
		x$x <- x$x[lg]
		x$y <- x$y[lg]
		if (ncp) x$y.x0 <- x$y.x0[lg, ]
	}

	if (from.t0 & (!fromt0Call) & (Args$t0 < x$x[1])) {
		x$x <- c(Args$t0, x$x)
		x$y <- c(0L, x$y)
		if (ncp) x$y.x0 <- rbind(0L, x$y.x0)
	}

	if ((!to.T) & toTCall & (x$x[length(x$x)] > t2)) {
		lg <- (x$x <= t2)
		x$x <- x$x[lg]
		x$y <- x$y[lg]
		if (ncp) x$y.x0 <- x$y.x0[lg, ]
	}

    	ymax <- max(x$y)

    	h <- pin.height - dp.height - gp.h - ti.h
    	A <- - ti.h * ymax/h
    	B <- (1 + (dp.height + gp.h)/h) * ymax
    	y2 <- ((B-A) + 1.08*(A+B))/2.16
    	y1 <- A + B - y2

      plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, ylim = c(y1, y2), axes = FALSE)
    	box()    

    	y3 <- par("usr")[4]
    	if (dp.legend) {
		dp.w <- sapply(dp.labels, strwidth, cex = dp.legend.cex)
	      dp.h <- sapply(dp.labels, strheight, cex = dp.legend.cex)	
	  	if (logic){
		 	text(par("usr")[1], par("usr")[4] - 0.5*max(dp.h[1:2]) - 0.1*par("cxy")[2], parse(text = paste("list(", 
				as.character(dp.labels[[1]]), ", ~~", as.character(dp.labels[[2]]), ")", sep="")), adj = 0, cex = dp.legend.cex)
	  	 	text(par("usr")[1], par("usr")[4] - max(dp.h[1:2]) - 0.5*dp.h[3] - 0.2*par("cxy")[2], dp.labels[[3]], adj = 0, cex = dp.legend.cex)
		 	y3 <- y3 - max(dp.h[1:2]) - dp.h[3] - 0.3*par("cxy")[2]		 
	  	}
	  	else{
		 	text(par("usr")[1], par("usr")[4] - 0.5 * dp.h[1] - 0.1*par("cxy")[2], dp.labels[[1]], adj = 0, cex = dp.legend.cex)
		 	text(par("usr")[1], par("usr")[4] - dp.h[1] - 0.5*dp.h[2] - 0.2*par("cxy")[2], parse(text = paste("paste(phantom(\" Diffusion process: \"), ~~", 
				as.character(dp.labels[[2]]), ")", sep="")), adj = 0, cex = dp.legend.cex)
		 	text(par("usr")[1], par("usr")[4] - sum(dp.h[1:2]) - 0.5*dp.h[3] - 0.3*par("cxy")[2], dp.labels[[3]], adj = 0, cex = dp.legend.cex)
		 		y3 <- y3 - sum(dp.h[1:2]) - dp.h[3] - 0.4*par("cxy")[2]		 
	  	}
	  	abline(h = y3)          
    	}

    	ticks <- axTicks(2)
    	axis(2, at = ticks[ticks <= y3], las = 1, mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    	axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    	title(main = "Approximate First-Passage-Time Density Function Plot", line = 0.85)
    	if (instants) title(xlab = "t", line = 4) else title(xlab = "t", line = 1.6)
    	if (ylab) title(ylab = parse(text = "g[1](t)"), line = 3)
     
    	if (any(growth.points, instants)) {
        	Y <- attr(x, "summary.fptl")[[1]]$instants
        	if (growth.points) {
            	i <- which(Y[, 1] >= x$x[1])
            	if (length(i) > 0) {
                		segments(Y[i, 1], par("usr")[3], Y[i, 1], y3, col = "darkgray", lwd = 1)
                		text(Y[i, 1] - 0.1 * par("cxy")[1], y3 - 0.5*strheight(expression(t), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
				parse(text = paste("t[~ ", i, "]", sep = "")), adj = 1, col = "darkgray")
            	}
      	}
		if (instants) {
            	x.t <- matrix(Y[, c(2, 3, 5)], ncol = 3)
            	segments(x.t, par("usr")[3], x.t, y3, lty = 8, lwd = 1)
            	text(Y[, 2] - 0.1 * par("cxy")[1], - 0.65*strheight(expression(t[1]), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2],
				parse(text = paste("t[~ ", 1:nrow(x.t), "]^{~ ", shQuote("*"), "}", sep = "")), adj = 1)
            	mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, line = 1.6, at = Y[, 3], 
				adj = 0, ...)
            	mtext(parse(text = paste("t[list(~ max,", 1:nrow(Y), ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, line = 3, at = Y[, 5], 
				adj = 0, ...)
        	}
    	}

	if (ncp){ 
		ymax <- max(x$y.x0)

    		h <- pin.height - dp.height - gp.h - ti.h
    		A <- - ti.h * ymax/h
    		B <- (1 + (dp.height + gp.h)/h) * ymax
    		y2 <- ((B-A) + 1.08*(A+B))/2.16
    		y1 <- A + B - y2
		dev.new()
		par(cex = 1, ps = 9)
		par(mar = c(3 + 2 * instants, 3.25 + ylab, 2.5, 1) + 0.1, ...)
    		matplot(x$x, x$y.x0, xlab = "", ylab = "", type = "l", lty = 1, las = 1, ylim = c(y1, y2), axes = FALSE)
   		box()    

    		y3 <- par("usr")[4]
    		if (dp.legend) {
			dp.w <- sapply(dp.labels, strwidth, cex = dp.legend.cex)
		      dp.h <- sapply(dp.labels, strheight, cex = dp.legend.cex)	
	  		if (logic){
		 		text(par("usr")[1], par("usr")[4] - 0.5*max(dp.h[1:2]) - 0.1*par("cxy")[2], parse(text = paste("list(", 
					as.character(dp.labels[[1]]), ", ~~", as.character(dp.labels[[2]]), ")", sep="")), adj = 0, cex = dp.legend.cex)
	  	 		text(par("usr")[1], par("usr")[4] - max(dp.h[1:2]) - 0.5*dp.h[3] - 0.2*par("cxy")[2], dp.labels[[3]], adj = 0, cex = dp.legend.cex)
		 		y3 <- y3 - max(dp.h[1:2]) - dp.h[3] - 0.3*par("cxy")[2]		 
	  		}
	  		else{
		 		text(par("usr")[1], par("usr")[4] - 0.5 * dp.h[1] - 0.1*par("cxy")[2], dp.labels[[1]], adj = 0, cex = dp.legend.cex)
		 		text(par("usr")[1], par("usr")[4] - dp.h[1] - 0.5*dp.h[2] - 0.2*par("cxy")[2], parse(text = paste("paste(phantom(\" Diffusion process: \"), ~~", 
					as.character(dp.labels[[2]]), ")", sep="")), adj = 0, cex = dp.legend.cex)
		 		text(par("usr")[1], par("usr")[4] - sum(dp.h[1:2]) - 0.5*dp.h[3] - 0.3*par("cxy")[2], dp.labels[[3]], adj = 0, cex = dp.legend.cex)
		 		y3 <- y3 - sum(dp.h[1:2]) - dp.h[3] - 0.4*par("cxy")[2]		 
	  		}
		}

	  	abline(h = y3)

		ticks <- axTicks(2)
    		axis(2, at = ticks[ticks <= y3], las = 1, mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    		axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    		title(main = "Approximate Conditional First-Passage-Time Density Functions Plot", line = 0.85)
		#  for selected values of the initial distribution
		title(xlab = "t", line = 1.6)
		if (ylab) title(ylab = parse(text = "g[1](t)"), line = 3)
	}
}
