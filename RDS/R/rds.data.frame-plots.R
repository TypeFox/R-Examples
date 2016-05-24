#' Diagnostic plots for the RDS recruitment process
#' @param x An rds.data.frame object.
#' @param plot.type the type of diagnostic.
#' @param stratify.by A factor used to color or stratify the plot elements.
#' @param ... Additional arguments for the underlying plot function if applicable.
#' @export
#' @details Several types of diagnostics are supported by the plot.type argument.
#' 'Recruitment tree' displays a network plot of the RDS recruitment process.
#' 'Network size by wave' monitors systematic changes is network size based on how far subjects are from the seed
#' 'Recruits by wave' displays counts of subjects based on how far they rare from their seed.
#' 'Recruit per seed' shows the total tree size for each seed.
#' 'Recruits per subject' shows counts of how many subjects are recruited by each subject who are non-terminal.
#' @return Either nothing (for the recruitment tree plot), or a ggplot2 object.
#' @method plot rds.data.frame
#' @examples 
#' data(fauxmadrona)
#' \dontrun{
#' plot(fauxmadrona)
#' }
#' plot(fauxmadrona, plot.type='Recruits by wave')
#' plot(fauxmadrona, plot.type='Recruits per seed')
#' plot(fauxmadrona, plot.type='Recruits per subject')
#' 
#' plot(fauxmadrona, plot.type='Recruits by wave', stratify.by='disease')
#' plot(fauxmadrona, plot.type='Recruits per seed', stratify.by='disease')
#' plot(fauxmadrona, plot.type='Recruits per subject', stratify.by='disease')
plot.rds.data.frame <- function(x,
		plot.type=c("Recruitment tree",
				"Network size by wave",					
				"Recruits by wave",
				"Recruits per seed",					
				"Recruits per subject"),
		stratify.by=NULL,
			...){
	
	#for R CMD check
	wave <- network <- ..n.. <- ..y.. <- seed <- color <- NULL
	ggplotGrob <- function(x) {ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))}
	
	x <- as.rds.data.frame(x)
	plot.type <- plot.type[1]
	tmp <- data.frame(id = x[[attr(x,"id")]],
			recruiter.id = x[[attr(x,"recruiter.id")]],
			network = x[[attr(x,"network.size.variable")]],
			wave = factor(get.wave(x)),
			seed = get.seed.id(x),
			color = if(is.null(stratify.by)) 0 else x[[stratify.by]]
	)
	p <- invisible()
	if(is.null(stratify.by)){
		if(plot.type=="Recruitment tree"){
			reingold.tilford.plot(x, 
					stratify.by=stratify.by,...)
			return(invisible())
		}else if(plot.type=="Network size by wave"){
			p <- ggplot(data=tmp) +
							geom_point(aes(x = wave,y = network,size = ..n..),stat = 'sum') +
							scale_size_area(name = '# of Subjects\nwith Identical \nValues')+
							stat_summary(aes(x = wave,y = network,
											ymax=..y..,ymin=..y..),fun.data = mean_sdl,
									geom = 'crossbar',color="red") +
							labs(title=plot.type)
		}else if(plot.type=="Recruits by wave"){
			p <- ggplot() +
							geom_bar(aes(x = as.factor(wave)),data=tmp) + xlab("wave") +
							labs(title=plot.type)
		}else if(plot.type=="Recruits per seed"){
			p <- ggplot() + geom_bar(aes(x = as.factor(seed)),data=tmp) + xlab("seed") + 
							labs(title=plot.type)
		}else if(plot.type=="Recruits per subject"){
			sd <- tmp$recruiter.id == get.seed.rid(x)
			p <- qplot(x=factor(as.vector(table(tmp$recruiter.id[!sd]))[-1]),
									xlab="# of Recruits") + 
							labs(title=plot.type)
		}
	}else{
		if(plot.type=="Recruitment tree"){
			reingold.tilford.plot(x, 
					vertex.color=stratify.by,...)
			return(invisible())
		}else if(plot.type=="Network size by wave"){
			p <- ggplot(data=tmp) +
							geom_point(aes(x = wave,y = network,colour = as.factor(color),size = ..n..),
									alpha = 0.5,stat = 'sum') +
							scale_size_area(name = '# of Subjects\nwith Identical \nValues')+
							stat_summary(aes(x = wave,y = network,colour = as.factor(color),ymax=..y..,ymin=..y..),
									fun.data = mean_sdl,geom = 'crossbar')+
							scale_colour_hue(name=stratify.by) + 
							labs(title=plot.type)
		}else if(plot.type=="Recruits by wave"){
			p <- ggplot2::ggplot() +
							geom_bar(aes(x = wave,fill=as.factor(color)),data=tmp)+
							scale_fill_hue(name=stratify.by) + 
							labs(title=plot.type)
		}else if(plot.type=="Recruits per seed"){
			p <- ggplot() + 
							geom_bar(aes(x = as.factor(seed),fill=as.factor(color)),data=tmp)+
							scale_fill_hue(name=stratify.by)+xlab("seed") +
							labs(title=plot.type)
		}else if(plot.type=="Recruits per subject"){
			rids <- get.rid(x)
			sd <- rids == get.seed.rid(x)
			tab <- table(rids[!sd])
			rids <- names(tab)
			idmap <- match(rids,get.id(x))
			col <- tmp$color[idmap]
			dat <- data.frame(Var2=col,value=as.vector(tab))
			p <- qplot(x=factor(dat$value),fill=factor(dat$Var2), xlab="# of Recruits") + 
							scale_fill_hue(guide=guide_legend(title=stratify.by)) + 
							labs(title=plot.type)
		}		
	}
	p
}

#' Plots the recruitment network using the Reingold Tilford algorithm.
#' @param x An rds.data.frame
#' @param vertex.color The name of the categorical variable in x to color the points with.
#' @param vertex.color.scale The scale to create the color palette.
#' @param vertex.size The size of the vertex points. either a number or the name of a 
#' column of x.
#' @param vertex.size.range If vertex.size represents a variable, vertex.size.range is a 
#' vector of length 2 representing the minimum and maximum cex for the points.
#' @param edge.arrow.size The size of the arrow from recruiter to recruitee.
#' @param vertex.label The name of a variable to use as vertex labels. NA implies no labels.
#' @param vertex.label.cex The size expansion factor for the vertex.labels.
#' @param vertex.frame.color the color of the outside of the vertex.points.
#' @param show.legend If true and either vertex.color or vertex.size represent variables, 
#' legends will be displayed at the bottom of the plot.
#' @param plot Logical, if TRUE then a plot is produced of recruitment tree.
#' ratio statistic with the observed statistics plotted as a vertical dashed line.
#' @param ... Additional parameters passed to plot.igraph.
#' @return A two-column vector of the positions of the nodes in the recruitment tree.
#' @export
#' @examples 
#' \dontrun{
#' data(fauxmadrona)
#' data(faux)
#' reingold.tilford.plot(faux)
#' reingold.tilford.plot(fauxmadrona,vertex.color="disease")
#' }
reingold.tilford.plot <-function(x, 
		vertex.color=NULL,
		vertex.color.scale = hue_pal(),
		vertex.size=2,
		vertex.size.range=c(1,5),
		edge.arrow.size=0,
		vertex.label.cex=.2,
		vertex.frame.color=NA, 
		vertex.label = get.id(x),
		show.legend=TRUE,
		plot=TRUE,
		...){	
	x <- as.rds.data.frame(x)
	
	if(!is.null(vertex.color)){
			color.name <- vertex.color
			color.var <- factor(x[[vertex.color]])
			levs <- levels(color.var)
			ncol <- length(levs)
			cols <- vertex.color.scale(ncol)
			color <- cols[as.integer(color.var)]
	}else
		color <- NULL
	
	if(is.character(vertex.size)){
		vertex.size.name <- vertex.size
		vertex.size <- as.numeric(x[[vertex.size]])
	}else{
		vertex.size.name <- deparse(substitute(vertex.size))
	}
	if(length(vertex.size)>1){
		vrange <- range(vertex.size,na.rm=TRUE)
		vertex.size <-vertex.size + vrange[1]
		vertex.size <- vertex.size / vrange[2]
		vertex.size <- vertex.size*(vertex.size.range[2]-vertex.size.range[1]) + vertex.size.range[1]
		vertex.size[is.na(vertex.size)] <- 0
	}
	if(length(vertex.label)==1 && !is.na(vertex.label)){
		vertex.label <- as.character(x[[vertex.label]])
	}
	
	id <- get.id(x)
	rid <- get.rid(x)
	sid <- get.seed.rid(x)
	seeds <- get.seed.id(x)
	el <- cbind(rid,id,seeds)
	el <- el[rid!=sid,]
	el <- el[order(el[,3]),, drop=FALSE]
	
	########
	#	generate layouts for each subgraph
	########	
	xyl <-list()
	grl <- list()
	for(seed in unique(seeds)){
		els <- el[el[,3]==seed, , drop=FALSE]
		if(nrow(els)>0){
			gr <-igraph::graph.edgelist(els[,1:2, drop=FALSE])
			lo <- igraph::layout.reingold.tilford(gr,root=seed)#,circular=TRUE)
			#sc <- mean(diff(unique(round(sort(lo[,1]),4))))
			tmp <- lo
			tmp[,1] <- round(lo[,1]/.25)
			overplt <- duplicated(tmp)
			if(any(overplt)){
				lo[overplt,2] <- lo[overplt,2] - .5
				#lo[overplt,1] <- lo[overplt,1] - .2
			}
		}else if(nrow(els)==1){
			gr <-igraph::graph.edgelist(els[,1:2, drop=FALSE])
			lo <- matrix(c(1,0,0,1),ncol=2)
		}else{
			gr <- igraph::graph.empty() + seed
			lo <- matrix(c(0,0),ncol=2)
		}
		if(!is.null(color)){
			i <- match(igraph::V(gr)$name,id)
			igraph::V(gr)$color <- color[i]
		}
		if(length(vertex.size)>1){
			i <- match(igraph::V(gr)$name,id)
			igraph::V(gr)$size <- vertex.size[i]
		}
		if(length(vertex.label)>1){
			i <- match(igraph::V(gr)$name,id)
			igraph::V(gr)$label <- vertex.label[i]			
		}
		
		xyl[[seed]] <- lo
		grl[[seed]] <- gr
	}


	########
	#	now layout subgraphs together without overlapping
	#	see: similar to wordcloud package
	########
	last <- 1
	overlap <- function(x1, y1, sw1, sh1, boxes) {
		s <- 0
		if (length(boxes) == 0) 
			return(FALSE)
		for (i in c(last,1:length(boxes))) {
			bnds <- boxes[[i]]
			x2 <- bnds[1]
			y2 <- bnds[2]
			sw2 <- bnds[3]
			sh2 <- bnds[4]
			if (x1 < x2) 
				overlap <- x1 + sw1 > x2-s
			else 
				overlap <- x2 + sw2 > x1-s
			
			if (y1 < y2) 
				overlap <- overlap && (y1 + sh1 > y2-s)
			else 
				overlap <- overlap && (y2 + sh2 > y1-s)
			if(overlap){
				last <<- i
				return(TRUE)
			}
		}
		last <<- 1
		FALSE
	}
	ord <- order(sapply(xyl,nrow),decreasing=TRUE)
	tstep=.1 
	rstep=.1
	
	boxes <- list()
	xyl2 <- xyl
	for(ind in 1:length(xyl)){
		i <- ord[ind]
		r <-0
		theta <- stats::runif(1,0,2*pi)
		x1 <- xo <- 0
		y1 <- yo <- 0
		wid <- diff(range(xyl[[i]][,1], na.rm=TRUE)) * 1.1
		sdx <- 1
		ht <- diff(range(xyl[[i]][,2], na.rm=TRUE)) * 1.1
		sdy <- 1
		isOverlaped <- TRUE
		while(isOverlaped){
			if(!overlap(x1-.5*wid,y1-.5*ht,wid,ht,boxes)){
				boxes[[length(boxes)+1]] <- c(x1-.5*wid,y1-.5*ht,wid,ht)
				isOverlaped <- FALSE
			}else{
				theta <- theta+tstep
				r <- r + rstep*tstep/(2*pi)
				x1 <- xo+sdx*r*cos(theta)
				y1 <- yo+sdy*r*sin(theta)
			}
		}
		xyl2[[i]][,1] <- xyl2[[i]][,1] + x1 - wid*.45 - min(xyl2[[i]][,1])
		xyl2[[i]][,2] <- xyl2[[i]][,2] + y1 - ht*.45 - min(xyl2[[i]][,2])
	}
	
	###########
	# now plot it
	###########
	t <- do.call(rbind,xyl2)
	if(plot){
	igr <- igraph::graph.disjoint.union(grl)
	nm <- do.call(c,lapply(grl,function(a)igraph::V(a)$name))
	vcol <- do.call(c,lapply(grl,function(a)igraph::V(a)$color))
	if(length(vertex.size)>1)
		vsize <- do.call(c,lapply(grl,function(a)igraph::V(a)$size))
	else
		vsize <- vertex.size
	if(length(vertex.label)>1){
		vlab <- do.call(c,lapply(grl,function(a)igraph::V(a)$label))
	}else
		vlab <- vertex.label
	igraph::V(igr)$name <- nm
	if(is.null(color)){
		igraph::plot.igraph(igr,
			layout=t,
			vertex.size=vsize,
			edge.arrow.size=edge.arrow.size,
			vertex.label.cex=vertex.label.cex,
			vertex.frame.color=vertex.frame.color,
			vertex.label=vlab,
			...)
	}else{
		igraph::plot.igraph(igr,
				layout=t,
				vertex.size=vsize,
				edge.arrow.size=edge.arrow.size,
				vertex.label.cex=vertex.label.cex,
				vertex.frame.color=vertex.frame.color,
				vertex.color=vcol,
				vertex.label=vlab,
				...)	
		if(show.legend)
			graphics::legend("bottomleft",legend=levs,col=cols,pch=16,
				title=color.name,horiz=TRUE,box.col=NA)
	}
	if(length(vsize)>1 && show.legend){
		levs <- paste(vrange," ")
		s <- vertex.size.range
		lg <- graphics::legend("bottomright",legend=levs,pt.cex=c(0,0),pch=16,
				title=vertex.size.name,horiz=TRUE,box.col=NA, x.intersp = 2)		
		t <- lg$text
		s <- s/200
		graphics::symbols(x=t$x-.075,y=t$y,circles=s,add=TRUE,inches=FALSE,bg="SkyBlue2")
	}
	}
	invisible(t)
}
