nma.networkplot <- function(s.id, t.id, data, title = "", trtname,
	alphabetic = TRUE, weight.edge = TRUE, adjust.thick = 5, weight.node = TRUE, adjust.node.size = 10,
	node.col = "orange", edge.col = "black", text.cex = 1,
	adjust.figsizex = 1.1, adjust.figsizey = 1.1){
	if(missing(s.id) | missing(t.id)){
		stop("both study id and treatment id must be specified.")
	}
	if(!missing(data)){
		s.id <- eval(substitute(s.id), data, parent.frame())
		t.id <- eval(substitute(t.id), data, parent.frame())
	}

	unique.sid <- unique(s.id)
	nstudy <- length(unique.sid)
	sid.ori <- s.id
	for(s in 1:nstudy){
		s.id[sid.ori == unique.sid[s]] <- s
	}

	unique.tid <- sort(unique(t.id))
	ntrt <- length(unique.tid)
	if(ntrt <= 2) stop("there are less than 3 treatments, no need for network plot.")
	if(!missing(trtname) & alphabetic){
		trtname.order <- order(trtname)
		unique.tid <- unique.tid[trtname.order]
	}
	if(!missing(trtname) & !alphabetic){
		trtname.order <- 1:ntrt
	}
	if(missing(trtname)){
		trtname <- unique.tid
		trtname.order <- 1:ntrt
	}
	if(length(trtname) != ntrt){
		stop("the length of trtname do not match the data.")
	}
	## make treatment id to be 1 to ntrt
	tid.ori <- t.id
	for(t in 1:ntrt){
		t.id[tid.ori == unique.tid[t]] <- t
	}

	polar <- pi/2 - 2*pi/ntrt*(0:(ntrt - 1))
	x <- cos(polar)
	y <- sin(polar)
	plot(x, y, axes = FALSE, xlab="", ylab="", cex = 0.1,
		xlim = c(-adjust.figsizex, adjust.figsizex),
		ylim = c(-adjust.figsizey, adjust.figsizey),
		main = title)

	wt <- matrix(0, ntrt, ntrt)
	for(t1 in 2:ntrt){
		for(t2 in 1:(t1 - 1)){
			study.t1 <- s.id[t.id == t1]
			study.t2 <- s.id[t.id == t2]
			study.t1.t2 <- intersect(study.t1, study.t2)
			wt[t1, t2] <- length(study.t1.t2)
		}
	}
	wt <- c(wt)
	if(weight.edge == TRUE){
		wt.unique <- unique(wt[wt > 0])
		wtmin <- min(wt.unique)
		wtmax <- max(wt.unique)
		if(wtmin < wtmax){
			wt[wt > 0] <- round(1 + adjust.thick*(wt[wt > 0] - wtmin)/(wtmax - wtmin))
		}else{
			wt[wt > 0] <- 2
		}
	}else{
		wt[wt > 0] <- 2
	}
	wt <- matrix(wt, ntrt, ntrt)

	for(t1 in 2:ntrt){
		for(t2 in 1:(t1 - 1)){
			if(t1 != t2 & wt[t1, t2] > 0){
				lines(x = x[c(t1, t2)], y = y[c(t1, t2)],
					lwd = wt[t1, t2], col = edge.col)
			}
		}
	}
	if(weight.node){
		wt <- matrix(0, ntrt, ntrt)
		for(t1 in 2:ntrt){
			for(t2 in 1:(t1 - 1)){
				study.t1 <- s.id[t.id == t1]
				study.t2 <- s.id[t.id == t2]
				study.t1.t2 <- intersect(study.t1, study.t2)
				wt[t1, t2] <- length(study.t1.t2)
			}
		}
		wt <- wt + t(wt)
		wt <- colSums(wt)
		node.sizes <- 3 + (wt - min(wt))/(max(wt) - min(wt))*adjust.node.size
		points(x, y, pch = 20, cex = node.sizes, col = node.col)
	}else{
		points(x, y, pch = 20, cex = 3, col = node.col)
	}

	sides <- numeric(ntrt)
	eps <- 10^(-4)
	for(t in 1:ntrt){
		if((polar[t] <= pi/2 & polar[t] > pi/4) |
			(polar[t] < -5*pi/4 & polar[t] >= -3*pi/2)){
			sides[t] <- 3
		}
		if(polar[t] <= pi/4 & polar[t] >= -pi/4){
			sides[t] <- 4
		}
		if(polar[t] < -pi/4 & polar[t] > -3*pi/4){
			sides[t] <- 1
		}
		if(polar[t] <= -3*pi/4 & polar[t] >= -5*pi/4){
			sides[t] <- 2
		}
	}
	for(t in 1:ntrt){
		if(weight.node){
		text(x = x[t], y = y[t], labels = trtname[trtname.order[t]],
			cex = text.cex)
		}else{
		text(x = x[t], y = y[t], labels = trtname[trtname.order[t]],
			pos = sides[t], cex = text.cex)
		}
	}
}