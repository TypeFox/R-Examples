plot.fptl <-
function (x, sfptl, from.t0 = TRUE, to.T = TRUE, dp.legend = TRUE, 
    dp.legend.cex = 1, ylab = TRUE, growth.points = TRUE, instants = TRUE, 
    ...) 
{
    if (!is.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fptl")))

    if (missing(sfptl)){
	  if (growth.points | instants | !from.t0 | !to.T) {
		G <- growth.intervals(x$x, x$y)
        	if (is.null(G)) {
            	growth.points = FALSE
            	instants = FALSE
            	from.t0 = TRUE
            	to.T = TRUE
        	}
        	else sfptl <- summary(x)
		FPTLInstants <- sfptl[[1]]$instants
	  }
    }
    else{
    	  if (!is.summary.fptl(sfptl)) 
            stop(paste(sQuote("sfptl"), " object is not of class ", shQuote("summary.fptl")))
	
	  Args <- formals(summary.fptl)	  
	  a <- as.list(attr(sfptl, "Call")[[1]])[-1]
	  Args[names(a)] <- a   
	  Y <- summary(x, Args$zeroSlope, Args$p0.tol, Args$k)
        attr(Y, "Call") <- attr(sfptl, "Call")
        if (!identical(sfptl, Y)) 
            stop("the x and sfptl objects do not match up")
	  FPTLInstants <- sfptl[[1]]$instants
    }

    par(cex = 1, ps = 9)    
    par(mar = c(3 + 2 * instants, 3.25 + ylab, 2.5, 1) + 0.1, ...)    
    pin.height <- par("pin")[2] - sum(par("mai")[c(1,3)]*(par("cex")-1))

    if (dp.legend) {
	  Call <- attr(x, "Call")	  
        Args <- as.list(Call[c("t0","T","x0")])
	  S <- Call$S
	  env <- Call$env
	  if (!is.null(env)){
		if (is.call(env)) Args <- c(Args, eval(env))
		else Args <- c(Args, attr(x, "vars")[[as.character(env)]])	
	  }
        dp.labels <- vector("list", 3)
        dp.labels[[1]] <- parse(text = deparse(substitute(paste(" Diffusion process: ", list(group("{", list(X(t), paste(t, "  in  ", group("[", 
                list(t0, T), "]"))), "}"), ~~P(X(t0) == x0) == 1)), Args)))
	  logic <- sapply(Args, is.character)
	  Args[logic] <- sapply(Args[logic], function(x) parse(text = x)[[1]])
	  dpMean <- eval(parse(text = paste("substitute(", attr(x, "dp")$mean, ", Args[logic])", sep = "")))
	  dpVar <- eval(parse(text = paste("substitute(", attr(x, "dp")$var, ", Args[logic])", sep = "")))
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

    if (!from.t0) {
        if (x$x[1] < FPTLInstants[1, 1]) {
            lg <- (x$x >= FPTLInstants[1, 1])
            x$x <- x$x[lg]
            x$y <- x$y[lg]
        }
    }

    if (!to.T) {
        if (x$x[length(x$x)] > FPTLInstants[nrow(FPTLInstants), 5]) {
            j <- which(x$x >= FPTLInstants[nrow(FPTLInstants), 5])[1]
            j <- j + which.min(x$y[j:length(x$y)]) - 1
            x$x <- x$x[1:j]
            x$y <- x$y[1:j]
        }
    }

    h <- pin.height - dp.height - gp.h - ti.h
    A <- - ti.h/h
    B <- 1 + (dp.height + gp.h)/h
    y2 <- ((B-A) + 1.08*(A+B))/2.16
    y1 <- A + B - y2
    plot(x$x, x$y, xlab = "", ylab = "", type = "l", las = 1, ylim = c(y1, y2), axes = FALSE)
    box()
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), las = 1, mgp = c(3, 0.35 + 0.4/sqrt(par("cex")), 0), tcl = -0.35)
    axis(1, mgp = c(3, 0.35, 0), tcl = -0.35)
    title(main = "First-Passage-Time Location Function Plot", line = 0.85)
    if (instants) title(xlab = "t", line = 4) else title(xlab = "t", line = 1.6)
    if (ylab) title(ylab = "FPTL(t)", line = 3)

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
    
    if (growth.points) {		
         segments(FPTLInstants[, 1], par("usr")[3], FPTLInstants[, 1], y3, col = "darkgray", lwd = 1)
	   text(FPTLInstants[, 1] - 0.1 * par("cxy")[1], y3 - 0.5*strheight(expression(t), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
			parse(text = paste("t[~ ", 1:nrow(FPTLInstants), "]", sep = "")), adj = 1, col = "darkgray")
    }

    if (instants) {
         x.t <- matrix(FPTLInstants[, c(2, 3, 5)], ncol = 3)
         y <- matrix(sfptl[[1]]$FPTLValues[, c(2, 3, 5)], ncol = 3)
         points(x.t, y, type = "p", cex = min(1, 1/par("cex")))
         segments(x.t, par("usr")[3], x.t, y, lty = 8, lwd = 1)
         segments(par("usr")[1], y, x.t, y, lty = 8, lwd = 1)
	   text(FPTLInstants[, 2] - 0.1 * par("cxy")[1], - 0.65*strheight(expression(t[1]), cex = sqrt(par("cex"))) - 0.1*par("cxy")[2], 
			parse(text = paste("t[~ ", 1:nrow(x.t), "]^{~ ", shQuote("*"), "}", sep = "")), adj = 1)
	   mtext(parse(text = paste("t[list(~ max,", 1:nrow(FPTLInstants), ")]^{~ ", shQuote("-"), "}", sep = "")), side = 1, line = 1.6, at = FPTLInstants[, 3], 
			adj = 0, ...)
         mtext(parse(text = paste("t[list(~ max,", 1:nrow(FPTLInstants), ")]^{~ ", shQuote("+"), "}", sep = "")), side = 1, line = 3, at = FPTLInstants[, 5], 
			adj = 0, ...)
    }    
}
