#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## the generic with x, y, ... comes from graphics
setMethod("plot", signature(x = "EMM", y = "missing"),
        function(x, y, method = c("igraph", "interactive", 
			"graph", "MDS", "cluster_counts",
                        "transition_counts"), data = NULL, 
                parameter=NULL, ...){ 

            method <- match.arg(method)

            p <- .get_parameters(list(
                            cluster_counts=TRUE,
                            arrow_width=TRUE,
                            arrows = "counts",      ## or "probabilities"
                            arrow_width_multiplier=1,
                            state_size_multiplier=1,
                            add_labels = TRUE,
                            cluster_labels = NULL,
                            
			    ## MDS
			    mark_clusters = TRUE,
                            draw_threshold = FALSE,
                            
			    ## graph, igraph
                            mark_states = NULL,
			    ## may be a vector of same length as mark_state
			    mark_states_color = "red", 
			    mark_transitions = NULL,
			    mark_transitions_color = "red",
			    
			    
			    nAttrs = list(),
                            eAttrs = list()
                            ), parameter)

	    if(size(x)<1) {
		warning("Empty EMM. No plot produced!")
		return(invisible(NULL))
	    }
                
	    ## check if we can use graphviz
	    if(method=="graph" && !.installed("Rgraphviz")) {
		warning("Rgraphviz not available/installed. Reverting to igraph!", 
			call. = FALSE)
		method <- "igraph"
	    }

            emm_centers <- cluster_centers(x)

            if(method=="cluster_counts") {
                pl <- barplot(sort(cluster_counts(x), decreasing=TRUE), 
                        ylab="Count", xlab = "State", ...)

            }else if(method=="transition_counts") {
                tr <- transitions(x)
                cnt <- transition(x, tr, type="counts")
                names(cnt) <- apply(tr, MARGIN=1, FUN = 
                        function(x) paste(x, collapse="->"))
                pl <- barplot(sort(cnt, decreasing=TRUE), 
                        ylab="Transition counts", ...)

            }else if(method=="igraph" || method=="interactive") {
                g <- as.igraph(x)
		
		if(method=="interactive") plot_fun <- tkplot
		else plot_fun <- plot

		if(p$arrow_width) {
		    e.width <- map(get.edge.attribute(g, "weight"),
			    c(.5,3))*p$arrow_width_multiplier
		}else{
		    e.width <- 1*p$arrow_width_multiplier
		}

		if(p$cluster_counts) {
		    v.size <- map(cluster_counts(x), 
			    c(8,30)) * p$state_size_multiplier
		}else{
		    v.size <- 10*p$state_size_multiplier
		}
		#v.cex <- v.size/10

		## cluster labels
		v.labels <- states(x)
		if(!is.null(p$cluster_labels)) 
		    v.labels <- p$cluster_labels
		if(!p$add_labels) 
		    v.labels <- ""

	
		## mark_states
		v.label.color <- "black"
                V(g)$color <- "white"
                if(!is.null(p$mark_states)) {
		    V(g)[p$mark_states]$color <- p$mark_states_color
                }

		## mark_transitions
		E(g)$color <- "black"
		if(!is.null(p$mark_transitions)) {
		E(g)$color <- "grey"
		    if(length(p$mark_transitions_color)==1) 
			p$mark_transitions_color <- rep(
			    p$mark_transitions_color, 
			    length(p$mark_transitions))
		    for(i in 1:length(p$mark_transitions)) {
			x <- strsplit(p$mark_transitions[i], "->")[[1]]
			E(g)[V(g)[x[1]] %->% V(g)[x[2]]]$color <- p$mark_transitions_color[i]
		    }
		}



		pl <- plot_fun(g, 
			...,
			
			## our default values
			layout=layout.fruchterman.reingold,
			#layout=layout.reingold.tilford(g, root=0)*-1,
			vertex.label.family="Helvetica",
			edge.label.family="Helvetica",
			#vertex.shape=v.shape,
			vertex.label=v.labels,
			vertex.size=v.size,
			
			## this does not work for interactive
			#vertex.label.cex=v.cex,
			
			## b/w
			#vertex.label.color="black",
			#vertex.color = "white",
			#edge.color = "black",
			vertex.label.color=v.label.color,
			#vertex.color = v.color,
			
			#edge.color = e.color,
			edge.width=e.width,
			#edge.label=e.label,
			#edge.label.cex=p$cex*.6,
			
			
			## only accepts a single value for now!
			#edge.arrow.size=(e.width-.5)*.3,
			edge.arrow.size=.5,
			asp=0 ## fill whole plot
			)


            }else if(method=="graph") {
		## this is Rgraphviz

                g <- as.graph(x)
		
		nAttrs <- p$nAttrs
                eAttrs <- p$eAttrs

		## vertex colors
		if(!is.null(p$mark_states)) {
                    nAttrs$fillcolor <- rep(p$mark_states_color, 
			    length(p$mark_states))
                    names(nAttrs$fillcolor) <- p$mark_states
                }

                ## vertex labels
                if(!is.null(p$cluster_labels)) {
                    names(p$cluster_labels) <- states(x)
                    nAttrs$label <- p$cluster_labels
                }

                if(!p$add_labels) {
                    p$cluster_labels <- rep("", size(x))
                    names(p$cluster_labels) <- states(x)
                    nAttrs$label <- p$cluster_labels
                }

                ## vertex size
                if(p$cluster_counts){
                    nAttrs$width <- (.5 +
                    cluster_counts(x)/max(cluster_counts(x))) * p$state_size_multiplier * .75
                }

		if(p$arrow_width && graph::numEdges(g)>0) {
		    ## this should work but the order in graph
		    ## lwd ordering seems to be broken in graph
		    #edges <- transitions(x)
		    
		    edg <- graph::edges(g)
		    from <- character()
		    to <- unlist(edg)
		    for(n in names(edg)) 
			from <- c(from, rep(n, length(edg[[n]])))
		    edges <- cbind(from=from, to=to)
		    ## end work around for graph

		    lwd <- transition(x, edges, type=p$arrows)
		    
		    ## normalize 
		    lwd <- lwd-min(lwd)
		    if(max(lwd)>0) lwd <- lwd/max(lwd)
		    lwd <- (1 + lwd * 4 ) * p$arrow_width_multiplier

		    names(lwd) <- apply(edges, 
		    	    MARGIN=1, FUN = function(z) paste(z, collapse="~"))
		    eAttrs$lwd <- lwd
		}

		## edge color
		if(!is.null(p$mark_transitions)) {
		    e.cols <- p$mark_transitions_color

		    if(length(e.cols)==1)
			e.cols <- rep(e.cols, length(p$mark_transitions))
		    names(e.cols) <- sapply(strsplit(p$mark_transitions, 
				    "->"), paste, collapse="~")
		
		    eAttrs$color <- e.cols
		}

                pl <- Rgraphviz::plot(g, recipEdges="distinct",
                        nodeAttrs = nAttrs, edgeAttrs = eAttrs, ...)
            } else {
                if(nrow(emm_centers)<3) stop('Less than 3 centers! Use method="graph".')

                ## self transitions are not visible for these plots
                x <- remove_selftransitions(x)

                if(is.null(data)){

                    if(ncol(emm_centers)==2 && 
                            tolower(x@measure)=="euclidean" ) {
                        ## we need no MDS
                        mds <- list(points = emm_centers, GOF=c(0,1))
			pts <- mds$points
			rownames(pts) <- states(x)
			sub <- ''

		    }else{
                        d <- dist(emm_centers, method=x@distFun)
                        mds <- cmdscale(d, eig=TRUE, add=TRUE)
                        #mds <- isoMDS(d, trace=FALSE)
			pts <- mds$points
			dimnames(pts) <- list(states(x), 
				c("Dimension 1", "Dimension 2"))
			sub <- paste("These two dimensions explain",
				round(100 * mds$GOF[2], digits = 2),
				"% of the point variability.")
			#sub <- ""
		    }

		    ## start plotting
		    plot(pts, type="n", ..., sub=sub)

                    ## use cex for point size
                    cex <- 2
                    if(p$cluster_counts) cex <- 
                    (2+cluster_counts(x)/max(cluster_counts(x)) * 5 ) * p$state_size_multiplier

                    ## arrows
                    edges <- transitions(x)
		    arrows_fromto <- cbind(
                            from_x = pts[edges[,1],1], from_y = pts[edges[,1],2],
                            to_x = pts[edges[,2],1],to_y = pts[edges[,2],2]
                            )
                    
		    ## edge colors
		    edge_colors <- "grey"
		    if(!is.null(p$mark_transitions)) {
			edge_colors <- rep("grey", nrow(edges))
			if(length(p$mark_transitions_color)==1)
			    p$mark_transitions_color <- rep(
				p$mark_transitions_color,
				length(p$mark_transitions))
			
			sl <- strsplit(p$mark_transitions, "->")
			for(i in 1:length(sl)) {
			    edge_colors[apply(edges == sl[[i]], MARGIN=1, 
				    all)] <- p$mark_transitions_color[i]
			}
		    }


                    ## make arrows shorter so they do not cover the nodes 
                    nodeRad2 <- cbind(
                            x=(sapply(cex, FUN=function(cx)  
                                            strwidth("o", cex=cx/1.2))/2)^2,
                            y=(sapply(cex, FUN=function(cx)
                                            strheight("o", cex=cx/1.2))/2)^2)

                    x2 <- (arrows_fromto[,3]-arrows_fromto[,1])^2
                    y2 <- (arrows_fromto[,2]-arrows_fromto[,4])^2
                    z2 <- x2+y2
                    shorten <- cbind(
                            x1=sqrt(nodeRad2[edges[,1],1]/z2 * x2),
                            y1=sqrt(nodeRad2[edges[,1],2]/z2 * y2),
                            x2=sqrt(nodeRad2[edges[,2],1]/z2 * x2),
                            y2=sqrt(nodeRad2[edges[,2],2]/z2 * y2)
                            )

                    signs <- matrix(1, ncol=2, nrow=nrow(shorten))
                    signs[arrows_fromto[,1] > arrows_fromto[,3],1] <- -1
                    signs[arrows_fromto[,2] > arrows_fromto[,4],2] <- -1
                    signs <- cbind(signs, signs*-1)

                    shorten <- shorten * signs
                    shorten[is.nan(shorten)] <- 0 ## take care of too short arrows
                    arrows_fromto <- arrows_fromto + shorten

                    ## lwd for arrows
                    lwd <- 1
                    if(p$arrow_width) {
                        lwd <- transition(x, edges[,1], edges[,2], type=p$arrows)
                        ## normalize 
                        lwd <- lwd-min(lwd)
                        lwd <- lwd/max(lwd)
                        lwd <- (1 + lwd * 4) * p$arrow_width_multiplier
                    }

                    ## arrows whines about zero length arrows
                    suppressWarnings(
                            arrows(arrows_fromto[,1], arrows_fromto[,2], 
                                    arrows_fromto[,3],arrows_fromto[,4],
                                    length=0.15, col=edge_colors, angle=20, lwd=lwd)
                            )

                    ## overplot points and text
		    if(!is.null(p$mark_states)) {
			s <- states(x) %in% p$mark_states
			points(pts[s,,drop=FALSE], cex=cex[s], 
				col=p$mark_states_color, pch=19)
			#points(pts[!s,,drop=FALSE], cex=cex[!s])
		    } 
		    
		    points(pts, cex=cex,...)

                    if(p$add_labels) {
                        ## plot labels
                        if(is.null(p$cluster_labels)) labels <- states(x)
                        else labels <- p$cluster_labels
                            cex <- cex/1.7
                        ## make sure double digit labels fit
                        cex <- cex * (strwidth("8")/
                                apply(cbind(strwidth(labels), strheight(labels)), 
                                        MARGIN=1, max))
                        text(pts, labels=labels, cex=cex)
                    }
	
		    pl <- mds

                } else {
		    ### MDS + data
                    ## project state centers onto dataset
                    
		    ## remove NA rows (resets)
		    data <- data[!apply(data, FUN=function(x) all(is.na(x)), MARGIN=1),]
		    
		    if(ncol(emm_centers)==2 && 
                            tolower(x@measure)=="euclidean" ) {
                        ## we need no MDS
                        mds <- list(points = emm_centers, GOF=c(0,1))
                        centers <- emm_centers
                        allpoints <- data
			sub <- ''
                    
                    }else{

                        d <- dist(rbind(emm_centers, data), method=x@distFun)
                        mds <- cmdscale(d, eig=TRUE, add=TRUE)
                        centers <- mds$points[1:nrow(emm_centers),1:2]
                        allpoints <- mds$points[-c(1:nrow(emm_centers)),1:2]
			colnames(allpoints) <- c("Dimension 1", "Dimension 2")
			sub <- paste("These two dimensions explain",
				round(100 * mds$GOF[2], digits = 2),
				"% of the point variability.")
		    }
                    ## points
                    if(p$mark_clusters){
                        point_center <- find_clusters(x, data, 
				match_cluster="exact")
                        ## make state name integer for pch
                        pch <- as.integer(factor(point_center, 
                                        levels = states(x)))
                        pch_center <- 1:size(x)
                        
                        ## make sure we stay below 25 for pch
                        while(any(pch>25, na.rm=TRUE)) {
                            reduce <- pch>25
                            reduce[is.na(reduce)] <- FALSE
                            pch[reduce] <- pch[reduce] -25
                        }
                        while(any(pch_center>25)) pch_center[pch_center>25] <-
                            pch_center[pch_center>25] -25


                        plot(allpoints, col="grey", pch=pch, ...,
                                sub=sub)

                        ## plot points which do not belong to any state
                        if(any(is.na(point_center)))
                            points(allpoints[is.na(point_center),,drop=FALSE], 
                                col="blue", pch=1)

                        ## add ellipses
                        ## FIXME: does not work with d>2
			if(p$draw_threshold) {
			    for (i in 1:size(x)) {
				thr <- x@tnn_d$var_thresholds[i]
				loc <- cluster_centers(x)[i,]
				lines(ellipsePoints(thr, thr, loc=loc, n=90), 
					col = 1, lty=2)
			    }
			}


                    }else{
                        plot(allpoints, xlab="Dimension 1", ylab="Dimension 2", 
                                col="grey", ...,
                                sub= paste("These two dimensions explain",
                                        round(100 * mds$GOF[2], digits = 2),
                                        "% of the point variability."))
                    }

                    cex <- 1

                    ## use cex for point size (scale: 1...3)
                    if(p$cluster_counts) cex <- 1+cluster_counts(x)/max(cluster_counts(x))*2

                    ## centers
                    if(p$mark_clusters) points(centers, 
                            col="red", pch=pch_center, cex=cex)
                    else points(centers, col="red", pch="+", cex=cex)
                    #points(centers, col="red", pch=1:size(x), cex=cex)
                    if(p$add_labels) text(centers, labels=states(x), pos=3)
                }
		    pl <- mds

            }
        
	return(invisible(pl))
	}
        )



setMethod("plot", signature(x = "TRACDS", y = "missing"),
        function(x, y, ...){ 
                g <- as.igraph(x)
		plot(g, ...,
			layout=layout.fruchterman.reingold,
			#layout=layout.reingold.tilford(g, root=0)*-1,
			vertex.label.family="Helvetica",
			edge.label.family="Helvetica",
			#vertex.shape=v.shape,
			vertex.label=states(x),
			#vertex.size=v.size,

			## b/w
			vertex.label.color="black",
			vertex.color = "white",
			edge.color = "black",

			#edge.color = e.color,
			#edge.width=e.width,
			#edge.label=e.label,
			#edge.label.cex=p$cex*.6,
			)	
	})


setMethod("plot", signature(x = "tNN", y = "missing"),
	function(x, y, 
		method = c("direct", "MDS"),
		data=NULL, parameter=NULL, ...){ 

	    method <- match.arg(method)

	    p <- .get_parameters(list(
			    draw_threshold = FALSE
			    ), parameter)


	    if(nclusters(x)<1) {
		warning("No clusters. No plot produced!")
		return(invisible(NULL))
	    }

	    centers <- as.data.frame(cluster_centers(x))

	    if(method=="direct") {
		### plot for data.frame does plot or pairs (for dim >2)
		if(is.null(data)) {
		    plot(centers, ...)
		}else{
		    plot(rbind(as.data.frame(data), centers),
			    col=c(rep("gray", nrow(data)), 
				    rep(2, nrow(cluster_centers(x)))),
			    pch=c(rep(1, nrow(data)),
				    rep(3, nrow(cluster_centers(x)))), ...)
				}


		## add ellipses?
		if(p$draw_threshold) {
		    for (i in 1:nclusters(x)) {
			thr <- x@tnn_d$var_thresholds[i]
			loc <- cluster_centers(x)[i,]
			lines(ellipsePoints(thr, thr, loc=loc, n=90), 
				col = 1, lty=3)
		    }
		}

	    }else if(method=="MDS") {

		if(is.null(data)) {
		    if(ncol(centers)==2 && 
			    tolower(x@measure)=="euclidean" ) {
			## we need no MDS
			mds <- list(points = centers, GOF=c(0,1))
			pts <- mds$points
			rownames(pts) <- clusters(x)
			sub <- ''

		    }else{
			d <- dist(centers, method=x@distFun)
			mds <- cmdscale(d, eig=TRUE, add=TRUE)
			#mds <- isoMDS(d, trace=FALSE)
			pts <- mds$points
			dimnames(pts) <- list(clusters(x), 
				c("Dimension 1", "Dimension 2"))
			sub <- paste("These two dimensions explain",
				round(100 * mds$GOF[2], digits = 2),
				"% of the point variability.")
		    }

		    plot(pts, ...,
			    sub= paste("These two dimensions explain",
				    round(100 * mds$GOF[2], digits = 2),
				    "% of the point variability."))

		} else {
		    ### MDS + data
		    ## project state centers onto dataset

		    ## remove NA rows (resets)
		    data <- data[!apply(data, FUN=function(x) 
			    all(is.na(x)), MARGIN=1),]

		    if(ncol(centers)==2 && 
			    tolower(x@measure)=="euclidean" ) {
			## we need no MDS
			mds <- list(points = centers, GOF=c(0,1))
			centerpoints <- centers
			allpoints <- data
			sub <- ''

		    }else{
			d <- dist(rbind(centers, data), method=x@distFun)
			mds <- cmdscale(d, eig=TRUE, add=TRUE)
			centerpoints <- mds$points[1:nrow(centers),1:2]
			allpoints <- mds$points[-c(1:nrow(centers)),1:2]
			colnames(allpoints) <- c("Dimension 1", "Dimension 2")
			sub <- paste("These two dimensions explain",
				round(100 * mds$GOF[2], digits = 2),
				"% of the point variability.")
		    }


		    plot(allpoints, col="grey", pch=1, ...,
			    sub=sub)

		 #   if(any(is.na(centerpoints)))
			points(centerpoints, col="red", pch=3)
		}

	    }else stop("Unknown method!")


	})
