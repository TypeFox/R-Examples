#ALG: plot.itree is derived from plot.rpart.

plot.itree <- function(x, uniform=FALSE, branch=1, compress=FALSE,
			     nspace, margin=0, minbranch=.3,highlight.color="black",do_node_re=FALSE, ...){
			     
    if(!inherits(x, "itree")) stop("Not an itree object")
    
	#ALG 6/12/2013:if we want to plot local risk estimates, then we need to do different
	#spacing between nodes. hence we add the flag 'do_node_re'.

	ext.methods <- c("regression_extremes","class_extremes")
	pur.methods <- 	c("regression_purity","class_purity")
	one.sided.methods <- c(ext.methods,pur.methods)

    if (!is.null(x$frame$splits)) x <- rpconvert(x)  #help for old objects
    if (nrow(x$frame) <= 1L)
        stop("fit is not a tree, just a root")

    #alg: 10/16/2012 if we're doing one-sided extremes, then enforce uniform=TRUE
    if(x$method %in% one.sided.methods){
    	# uniform explicitly specified FALSE is a problem.
		if(!missing(uniform) && uniform==FALSE){
			stop("cannot specify uniform=FALSE for one sided methods")
		}
		uniform <- TRUE
    }

	#check if a highlight.color is passed with a non one-sided method.
	#If so, stop and print error
	highlight.color <- tolower(highlight.color)
	if(highlight.color!="black" && !(x$method %in% one.sided.methods ))
	{
		stop("Non-black highlight.color is only valid for 'extremes' or 'purity' methods")
	}
	# check if the highlight.color is valid, if not stop and print an error.
	if( (!tolower(highlight.color) %in% colors())){

		stop(paste("Error in highlight.color : invalid color name",highlight.color))
	}


    if (compress & missing(nspace)) nspace <- branch
    if (!compress) nspace <- -1L     #means no compression

	#6/24. made to match current rpart version in which envir!=Global environment
	#similarly adjusted in itree.branch.
    assign(paste(".itree.parms", dev.cur(), sep = "."),
            list(uniform=uniform, branch=branch, nspace=nspace,
		 minbranch=minbranch), envir=itree_env)

    # ### define the plot region
    temp <- rpartco(x)
    xx <- temp$x
    yy <- temp$y

    #alg 6/13/2013. add the if() statement here because if we're going
    # to plot  local risk estimates, more space is needed.
    if(do_node_re==TRUE){
		minx <- min(xx); maxx <- max(xx)
		miny <- min(yy); maxy <- max(yy)
		xx2 <- xx+ (xx==minx)*(-.01) + (xx==maxx)*(.01)
		yy2 <- yy+ (yy==miny)*(-.1) + (yy==maxy)*(.01)
	    temp1 <- range(xx2) + diff(range(xx2))*c(-margin, margin)
	    temp2 <- range(yy2) + diff(range(yy2))*c(-margin, margin)
	    plot(temp1, temp2, type='n', axes=FALSE, xlab='', ylab='', ...)
	}
	else{
	# otherwise do what is in plot.rpart.s from the rpart package
	    temp1 <- range(xx) + diff(range(xx))*c(-margin, margin)
	    temp2 <- range(yy) + diff(range(yy))*c(-margin, margin)
	    plot(temp1, temp2, type='n', axes=FALSE, xlab='', ylab='', ...)
	}
	

    # Draw a series of horseshoes or V's, left son, up, down to right son
    #   NA's in the vector cause lines() to "lift the pen"
    node <- as.numeric(row.names(x$frame))
    temp <- itree.branch(xx, yy, node, branch)


    if (branch > 0) text(xx[1L], yy[1L], '|')

    # ALG 9/4/2012: do coloring for one-sided criteria
	if((x$method %in% one.sided.methods) && highlight.color!= "black")
	{
		d1 <- dim(temp$x)[1]
		d2 <-  dim(temp$x)[2]

		#split up so we can plot left and right branches separately
		temp$x <- rbind(temp$x[1:2,], (temp$x[2,]+temp$x[3,])/2,rep(NA,d2),
						(temp$x[2,]+temp$x[3,])/2,temp$x[3:d1,])
		temp$y <- rbind(temp$y[1:2,],temp$y[2,],rep(NA,d2),
						temp$y[3,],temp$y[3:d1,])

		#redefine dimension
		d1 <- dim(temp$x)[1]
		
		# alg: need to do the lines in different colors
		######### PURITY:
		if(x$method %in% pur.methods){
			#We print the purer side in highlight.color.
			#figure out which child node is purer on a per-observation basis...
			#left.purity <-  x$frame$dev[temp$leftnode.row] / x$frame$n[temp$leftnode.row]
			#right.purity <- x$frame$dev[temp$rightnode.row] / x$frame$n[temp$rightnode.row]
			#as of 3/27/2013 it's not per observation...
			left.purity <-  x$frame$dev[temp$leftnode.row]
			right.purity <- x$frame$dev[temp$rightnode.row]
			
			left.is.purer <- (left.purity <= right.purity)

			# When left branch is purer, left of horseshoe is highlight.color
			lines(c(temp$x[1:4,left.is.purer==TRUE]),
					c(temp$y[1:4,left.is.purer==TRUE]), col=highlight.color)
			# ... and the right is black.
			lines(c(temp$x[5:d1,left.is.purer==TRUE]),
					c(temp$y[5:d1,left.is.purer==TRUE]))

			# When right branch is purer, the left of the horseshoe is black
		    lines(c(temp$x[1:4,left.is.purer==FALSE]),
		    		c(temp$y[1:4,left.is.purer==FALSE]))
		    # ... and the right is highlight.color.
		    lines(c(temp$x[5:d1,left.is.purer==FALSE]),
		    		c(temp$y[5:d1,left.is.purer==FALSE]), col=highlight.color)
	    	invisible(list(x=xx, y=yy))
		}
		
		##### EXTREMES:
		if(x$method %in% ext.methods){
						
			##### EXTREMES + CLASSIFICATION
			if(class(x$parms)=="list"){ #extremes highlighting for classification
				coi <- x$parms$classOfInterest
				numlevel <- length(unique(x$y))
				#colnum <- 1+numlevel+coi
				colnum <- coi
				
				#the fraction of obs = classOfInterest for each node.  
				#left.branch.is.greater <- (x$frame[temp$leftnode.row,][,colnum] >= x$frame[temp$rightnode.row,][,colnum])

				#6/13/2014: changed the below to use the new naming of tree$frame
				# fraction of each class in each node: 
				yval2mat <- x$frame[, grep("wt.frac.class",colnames(x$frame))]
				#left.branch.is.greater <- (yval2mat[temp$leftnode.row,][,colnum] >= yval2mat[temp$rightnode.row,][,colnum])
				left.branch.is.greater <- (yval2mat[temp$leftnode.row,colnum] >= yval2mat[temp$rightnode.row,colnum])
								
				#high means: highlight branch generating the large mean in highlight.color

				# When left branch is greater, left of horseshoe is green
				lines(c(temp$x[1:4,left.branch.is.greater==TRUE]),
						c(temp$y[1:4,left.branch.is.greater==TRUE]), col=highlight.color)
				# ... and the right is black.
				lines(c(temp$x[5:d1,left.branch.is.greater==TRUE]),
						c(temp$y[5:d1,left.branch.is.greater==TRUE]))

				# When right branch is greater, the left of the horseshoe is black
			    lines(c(temp$x[1:4,left.branch.is.greater==FALSE]),
			    		c(temp$y[1:4,left.branch.is.greater==FALSE]))
			    # ... and the right is highlight.color.
			    lines(c(temp$x[5:d1,left.branch.is.greater==FALSE]),
			    		c(temp$y[5:d1,left.branch.is.greater==FALSE]), col=highlight.color)
		    	invisible(list(x=xx, y=yy))
											
			}
			###### EXTREMES + REGRESSSION
			else{  
				#figure out if the left branch has higher mean
				left.branch.is.greater <- (x$frame$yval[temp$leftnode.row] >= x$frame$yval[temp$rightnode.row])
						
				if(x$parms==1){
					#high means: highlight branch generating the large mean in highlight.color
	
					# When left branch is greater, left of horseshoe is green
					lines(c(temp$x[1:4,left.branch.is.greater==TRUE]),
							c(temp$y[1:4,left.branch.is.greater==TRUE]), col=highlight.color)
					# ... and the right is black.
					lines(c(temp$x[5:d1,left.branch.is.greater==TRUE]),
							c(temp$y[5:d1,left.branch.is.greater==TRUE]))
	
					# When right branch is greater, the left of the horseshoe is black
				    lines(c(temp$x[1:4,left.branch.is.greater==FALSE]),
				    		c(temp$y[1:4,left.branch.is.greater==FALSE]))
				    # ... and the right is highlight.color.
				    lines(c(temp$x[5:d1,left.branch.is.greater==FALSE]),
				    		c(temp$y[5:d1,left.branch.is.greater==FALSE]), col=highlight.color)
			    	invisible(list(x=xx, y=yy))
				}
				if(x$parms==-1){
					#low means: highlight branch generating the small mean  in highlight.color
	
					# When left branch is greater, left of horseshoe is black
					lines(c(temp$x[1:4,left.branch.is.greater==TRUE]),
							c(temp$y[1:4,left.branch.is.greater==TRUE]))
					# ... and the right is red.
					lines(c(temp$x[5:d1,left.branch.is.greater==TRUE]),
							c(temp$y[5:d1,left.branch.is.greater==TRUE]),col=highlight.color)
	
					# When right branch is greater, the left of the horseshoe is highlight.color
				    lines(c(temp$x[1:4,left.branch.is.greater==FALSE]),
				    		c(temp$y[1:4,left.branch.is.greater==FALSE]),col=highlight.color)
				    # ... and the right is black.
				    lines(c(temp$x[5:d1,left.branch.is.greater==FALSE]),
				    		c(temp$y[5:d1,left.branch.is.greater==FALSE]))
			    	invisible(list(x=xx, y=yy))
				}
			}
			#end of highlighting for extremes+regression
		}
	}
	#alg: otherwise do the usual/existing rpart.
	if(highlight.color=="black" || (!(x$method %in% one.sided.methods)) ){
	    lines(c(temp$x), c(temp$y))
	    invisible(list(x=xx, y=yy))
	    #return(list(x=temp$x, y=temp$y))
	}
}




