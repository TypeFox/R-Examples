# plot.r ##################################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .get.leaves				      Reads the labels of variables. (Internal function)
#  .X.coord.var				      Computes the x-coordinates for the variables. (Internal function)
#  .X.coord.par				      Computes the x-coordinates for the parameters. (Internal function)
#  .get.coord				        Supplementary function of .X.coord.par. (Internal function)
#  .already.coord			      Supplementary function of .X.coord.par. (Internal function)
#  .shift				 	          Computes the y-coordinates of the variables. (Internal function) 
#  .line					          Supplementary function of .plot.lines.circles. (Internal function)
#  .plot.lines.circles		  Plots the lines, circles and the names of the variables. (Internal function)
#  .plot.rectangles			    Plots the rectangles and included text. (Internal function)
#  .circle					        Supplementary function of .plot.lines.circles. (Internal function)
#  .rectangle				        Supplementary function of .plot.rectanges. (Internal function)
#  plot.hac					        Produces the plot of the HAC structure.
##########################################################################################################################

.get.leaves = function(tree){
	if(length(tree)==1){tree = tree[[1]]}
	rapply(tree, classes = "character", f = function(r)r, how = "unlist")
}

#-------------------------------------------------------------------------------------------------------------------------------

.X.coord.var = function(tree, s){
	if( (!("coord" %in% names(tree))) & (class(tree) != "character")){
		n = length(tree) - 1
			for(i in 1:n){
				tree[[i]] = .X.coord.var(tree[[i]], s = s)
	}}else{
		if(class(tree) == "character"){
			n = length(tree)
        	tree = list(tree, coord = list(x = integer(n), y = integer(n)))
        	for(i in 1:n){
        		tree$coord$x[i] = which((tree[[1]][i] == s))
        	}
    }}
	tree
}

#-------------------------------------------------------------------------------------------------------------------------------

.X.coord.par = function(tree){
	
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree)
	s = sapply(tree, .already.coord)[-n]
	
	if(any(s)){
		if(any(!s)){
			include = c(1:(n-1))[which(s)]
			exclude = c(1:(n-1))[which(!s)]
			tree[exclude] = tree.down = lapply(tree[exclude], FUN = .X.coord.par)
			if(length(tree.down)==1){tree.down = tree.down[[1]]}			
			m = c(mean(sapply(tree[include], .get.coord)), mean(sapply(tree.down, .get.coord)))
			tree[[n]] = as.list(tree[[n]])
			tree[[n]]$coord$x = mean(m)
			tree[[n]]$coord$y = 0
		}else{
			include = c(1:(n-1))[which(s)]
			m = sapply(tree[include], FUN = .get.coord)
			tree[[n]] = as.list(tree[[n]])
			tree[[n]]$coord$x = mean(m)
			tree[[n]]$coord$y = 0
	}}else{
		tree[-n] = lapply(tree[-n], FUN = .X.coord.par)
		if(length(tree[-n])==1){tree[-n] = tree[-n][[1]]}
        m = sapply(tree[-n], .get.coord)	
		tree[[n]] = as.list(tree[[n]])
		tree[[n]]$coord$x = mean(m)
		tree[[n]]$coord$y = 0
	}	
	tree
}

#-------------------------------------------------------------------------------------------------------------------------------

.get.coord = function(tree){
	if((class(tree[[1]]) == "character") | (class(tree[[1]]) == "numeric")){
		n = length(tree)
		mean(tree[[n]][[1]])
	}else{
	if(class(tree[[1]]) == "list"){
		mean(unlist(sapply(tree, FUN = .get.coord)))
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.already.coord = function(tree){("coord" %in% names(tree))}

#-------------------------------------------------------------------------------------------------------------------------------

.shift = function(tree){
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree)
	m = integer(1)
	s = integer(n-1)
	for(i in 1:(n-1)){s[i]=is.character(tree[[i]][[1]])}
	
	if(any(s == 0)){
		exclude = c(1:(n-1))[which(s==0)]
			for(i in 1:(n-1)){
				d = length(tree[[i]])
				if(is.character(tree[[i]][[1]])){
					tree[[i]][[2]][[2]] = rep(tree[[n]][[2]][[2]], length(tree[[i]][[2]][[2]])) - 1
				}else{
					tree[[i]][[d]][[2]][[2]] = array(tree[[n]][[2]][[2]], dim = length(tree[[i]][[d]][[2]][[2]])) - 1
			}}
			tree[exclude] = lapply(tree[exclude], .shift)
	}else{
		for(i in 1:(n-1)){
			tree[[i]][[2]][[2]] = rep(tree[[n]][[2]][[2]], length(tree[[i]][[2]][[2]])) - 1}
	}
	tree
}

#-------------------------------------------------------------------------------------------------------------------------------

.line = function(upper, lower, col, lwd, ...){
	lines(c(upper), c(lower), col = col, lwd = lwd, ...)
}	

#-------------------------------------------------------------------------------------------------------------------------------

.plot.lines.circles = function(tree, circles, bg, fg, col, col.t, lwd, ...){
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree)
	
	for(i in 1:(n-1)){
		if(class(tree[[i]][[1]])=="list"){
			.line(c(tree[[n]]$coord$x, tree[[i]][[length(tree[[i]])]]$coord$x), c(tree[[n]]$coord$y, tree[[i]][[length(tree[[i]])]]$coord$y), col = col, lwd = lwd, ...)
			.plot.lines.circles(tree[[i]], circles = circles, bg = bg, fg = fg, col = col, col.t = col.t, lwd = lwd, ...)
		}else{
			for(j in 1:length(tree[[i]]$coord$x)){
			.line(c(tree[[n]]$coord$x, tree[[i]]$coord$x[j]), c(tree[[n]]$coord$y, tree[[i]]$coord$y[j]), col = col, lwd = lwd, ...)
			.circle(a = tree[[i]]$coord$x[j], b = tree[[i]]$coord$y[j], L = tree[[i]][[1]][j], radius = circles, fg = fg, bg = bg, col = col, col.t = col.t, lwd = lwd, ...)
	}}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.plot.rectangles = function(tree.coord, tree, h = 0.45, l = 1.2, z, index, numbering, s.params, theta, type, digits, fg, bg, col, col.t, lwd = lwd, ...){
	if(length(tree)==1){tree = tree[[1]]}
    if(length(tree.coord)==1){tree.coord = tree.coord[[1]]}
	n = length(tree)
    nn = length(tree.coord)
    
    stopifnot(n==nn)
	
	for(i in 1:(n-1)){
		if(class(tree.coord[[i]][[1]])=="list"){
        n = length(tree[[i]])
			.rectangle(a = tree.coord[[i]][[n]]$coord$x, b = tree.coord[[i]][[n]]$coord$y, L = tree[[i]], l = l, h = h, z = z, index = index, numbering = numbering, s.params = s.params, theta = theta, type = type, digits = digits, fg = fg, bg = bg, col = col, col.t = col.t, lwd = lwd, ...)
			.plot.rectangles(tree.coord = tree.coord[[i]], tree = tree[[i]], h = h, l = l, z = z, index = index, numbering = numbering, s.params = s.params, theta = theta, type = type, digits = digits, fg = fg, bg = bg, col = col, col.t = col.t, lwd = lwd, ...)
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------
	
.circle = function(a, b, L, radius, fg, bg, col, col.t, lwd, ...){
	symbols(a, b, circles = radius, add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
	text(a, b, L, col = col.t, ...)
}

#-------------------------------------------------------------------------------------------------------------------------------

.rectangle = function(a, b, L, l, h, z, index, numbering, s.params, theta, type, digits, fg, bg, col, lwd, col.t, ...){
	n = length(L)
	if(theta){
		if(!index){
			symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			text(a, b, label = bquote(paste(theta == .(round(L[[n]][[1]], digits = digits)))), col = col.t, ...)
		}else{
            if(!numbering){
                symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			    text(a, b, label = bquote(paste(theta[.(.allocate.all(L, theta = FALSE))]) == .(round(L[[n]][[1]], digits = digits))), col = col.t, ...)
            }else{
                symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			    text(a, b, label = bquote(paste(theta[.(which(s.params == L[[n]][[1]]))]) == .(round(L[[n]][[1]], digits = digits))), col = col.t, ...)
    }}}else{
		if(!index){
			symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			text(a, b, label = bquote(paste(tau == .(round(theta2tau(L[[n]][[1]], type), digits = digits)))), col = col.t, ...)
		}else{
            if(!numbering){
                symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			    text(a, b, label = bquote(paste(tau[.(.allocate.all(L, theta = FALSE))]) == .(round(theta2tau(L[[n]][[1]], type), digits = digits))), col = col.t, ...)
            }else{
                symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, fg = fg, bg = bg, lwd = lwd, ...)
			    text(a, b, label = bquote(paste(tau[.(which(s.params == L[[n]][[1]]))]) == .(round(theta2tau(L[[n]][[1]], type), digits = digits))), col = col.t, ...)
    }}}
}

#--------------------------------------------------------------------------------------------------------------------------------

.min.y = function(tree){
	n = length(tree)
	if(n == 1){tree = tree[[1]]}
	
	n = length(tree)
	s = logical(n)
	m = integer(1)
	
	for(i in 1:n){
		if(class(tree[[i]][[1]])=="list"){
			s[i] = TRUE
		}else{
			s[i] = FALSE
	}}
	
	if(any(s)){
		 m = sapply(tree[which(s)], .min.y)
	}else{
		m = min(tree[which(!s)][[1]]$coord$y)
	}
	min(unlist(m))
}

#--------------------------------------------------------------------------------------------------------------------------------

plot.hac = function(x, xlim = NULL, ylim = NULL, xlab = "", ylab = "", col = "black", fg = "black", bg = "white", col.t = "black", lwd = 2, index = FALSE, numbering = FALSE, theta = TRUE, h = 0.4, l = 1.2, circles = 0.25, digits = 2, ...){

	tree = x$tree
	main.var = .get.leaves(tree)
	dd = length(main.var)
	s = 0.3 * dd
	if(is.null(xlim)){xlim = c(1 - circles, dd + circles)}
	s.params = get.params(x, sort.v = TRUE, decreasing = TRUE)

d = length(tree)
tree.coord = .shift(.X.coord.par(.X.coord.var(tree, main.var)))

if(is.null(ylim)){ylim = c((.min.y(tree.coord) - 2 * (circles + 0.05)), h / 2)}
	plot(x = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = FALSE, col = "white", ...)
	.plot.lines.circles(tree.coord, circles = circles, bg = bg, fg = fg, col = col, col.t = col.t, lwd = lwd, ...)
	.plot.rectangles(tree.coord = tree.coord, tree = tree, h = h, l = l, z = s, index = index, numbering = numbering, s.params = s.params, theta = theta, type = x$type, digits = digits, bg = bg, fg = fg, col = col, col.t = col.t, lwd = lwd, ...)
	.rectangle(a = tree.coord[[d]]$coord$x, b = tree.coord[[d]]$coord$y, L = tree, l = l, h = h, z = s, index = index, numbering = numbering, s.params = s.params, theta = theta, type = x$type, digits = digits, fg = fg, bg = bg, col = col, col.t = col.t, lwd = lwd, ...)
}
