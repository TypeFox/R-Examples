################################################################################
# PLOT.MANYLM, PLOT.MANYGLM:                                                   #
# Plot for evaluation of goodness of fit for lm.mvabund objects                #
################################################################################

default.plot.manylm  <- function(x, 
				which = 1:4, 
				caption = c("Residuals vs Fitted","Normal Q-Q", "Scale-Location", "Cook's distance"), 
				overlay=TRUE,
    				n.vars=Inf, 
				var.subset=NULL, 
				panel = if (add.smooth) panel.smooth else points,
    				sub.caption = NULL, 
				main = "", 
				ask, 
				...,
    				id.n = if (overlay) 0 else 3, 
				labels.id =  rownames(as.matrix(residuals(x))),
	    			cex.id = 0.75, 
				qqline = TRUE, 
				cook.levels = c(0.5, 1),
    				add.smooth = if(!is.null(getOption("add.smooth"))){ 
							getOption("add.smooth")
    						 } else TRUE, label.pos = c(4, 2), 
				cex.caption = 1, 
				asp = 1,
			    	legend.pos= if(length(col)==1) "none" else "nextplot",
    				mfrow= if(overlay) { 
						length(which)+(legend.pos=="nextplot")
    					 } else if(write.plot=="show") c(min(n.vars,3),length(which)) 
					 else length(which), 
				mfcol=NULL, write.plot="show", 
				filename="plot.mvabund",
    				keep.window= if(is.null(c(mfrow,mfcol))) TRUE 
						 else FALSE,  
				studentized=TRUE) 
{	


	allargs <- match.call(expand.dots = FALSE)
     	dots <- allargs$...
    
     	dev <- dev.list()
     	dev.name <- getOption("device")

#     	if(is.null(dev.name)) stop("Make sure that the 'device' option has a valid value,
#     		e.g. 'options(device = 'windows')'. Allowed values here are 'windows',
#			 'win.graph', 'x11', 'X11'.")

#     	if(!(any(dev.name == c("windows", "win.graph", "x11", "X11")) ))
#      	stop("Make sure that the 'device' option has a valid value,
#     			e.g. 'options(device = 'windows')'. Allowed values here are 
#				'windows', 'win.graph', 'x11', 'X11'.")
    
    	if (write.plot!="show") {
      	if (write.plot=="eps" | write.plot=="postscript") {
        		postscript(paste(filename,".eps", sep="") )
        	} else if (write.plot=="pdf") {
			pdf(paste(filename,".pdf", sep="") )
        	} else if (write.plot=="jpeg" ){
			jpeg(paste(filename,".jpeg", sep=""))
        	} else if (write.plot=="bmp" ){
			bmp(paste(filename,".bmp", sep=""))
        	} else if (write.plot=="png" ){
			png(paste(filename,".png", sep=""))
		}  
    		on.exit( dev.off() )
    	}

    	na.action       <- x$na.action
    	na.action.type  <- attr(na.action, "class")
    	if(!is.null(na.action))  message("Due to NA values the case(s) ", na.action, " were discarded.")

   	if (length(dots)>0) {  
		# Delete arguments in ... that are defined lateron and cannot be used twice
   		# in the plot function.
    		deactive <- c("xlab", "ylab", "ylim", "sub", "type") 
		deactivate <- (1:length(dots))[names(dots) %in% deactive ]
    
		for (i in length(deactivate):1) {
   	 		dots[ deactivate[i] ]<-NULL 				##fixed up [[]], due to compile error (v2.10).
		}

    		dots <- lapply( dots, eval, parent.frame() )
    		if( "col.main" %in% names(dots) ) colmain <- dots$col.main
    		else colmain <- par("col.main")

    	} else {
    		colmain <- par("col.main")
    	}  

    	if (!inherits(x, c("manylm", "manyglm"))) 
      	warning("use 'plot.manylm' only with \"manylm\" or \"manyglm\" objects")
        
    	if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        	stop("'which' must be in 1:4")
    
    	show <- rep.int(FALSE, times = 4)
    	show[which] <- TRUE
    	
	if(ncol(x$x) > 0 ) empty <- FALSE 
	else empty <- TRUE

    	if(empty && show[4]) 
      	if (length(which)==1) {
    		stop("Plot no. 4 cannot be drawn, as Cooks distance cannot be calculated for an empty model")
       	} else {
       		warning("Plot no. 4 cannot be drawn, as Cooks distance cannot be calculated for an empty model")
       		show[4] <- FALSE
       		which   <- which[which==4]
        }

    	r <- r.orig <- as.matrix(residuals(x))
    	yh <- predict(x)   # the linear predictor for glm's
    	
	w <- weights(x)

    	if(any(na.action.type == "exclude")) {
        	r <- r[- na.action, ]
        	r.orig <- r.orig[- na.action, ]
        	if(!is.null(w)) w <- w[- na.action]
       
        	yh <- yh[- na.action, ]
    	}
    

    	if(any(na.action.type == "pass") | is.null(na.action.type)) {
      	    which.na.pass <- which(is.na(r[1,]))      	
	    if(length(which.na.pass)>0){
       		yh  <-  yh[ , - which.na.pass] # !!! needs to be checked, maybe this line must be deleted !!!
       		r   <- r[ , - which.na.pass]   # the whole column with NA values in y will be NA in residuals
       		r.orig <- r.orig[, -  which.na.pass ]
            } 
	}
    
    	miss.varsubset <- missing(var.subset) | is.null(var.subset)
    	# Change logical var.subset to numerical var.subset, if necessary. Note that
    	# NA values are logical as well, but should be excluded here.
    
	if(!miss.varsubset){
    	    if(is.logical(var.subset) & any(!is.na(var.subset)))
                    var.subset <- which(var.subset[!is.na(var.subset)])
    	}
    
	miss.varsubset <- !is.numeric(var.subset) # If this function is called within
    								# another, the missing function could be tricked out.

	if(is.null(labels.id)) labels.id <- as.character(1:nrow(r))

    	if (!is.null(w)) {
         	wind <- w != 0
       		w  <- w[wind]
        	r <- r[wind,, drop=FALSE]
        	r.orig <- r.orig[wind,, drop=FALSE]
        	yh <- yh[wind,, drop=FALSE]  
        	labels.id <- labels.id[wind]
    	}
    
	n <- nrow(r)
    	p <- ncol(r)

    	########## BEGIN edit var.subset, n.vars and r & fitted values  ############
    	# subset allows double variables
    	var.subset.dim <- length(var.subset)
    	if (miss.varsubset){
           	if (n.vars>p) n.vars <- min(n.vars,p)

                y <- as.matrix(x$y)

       	        if(any(na.action.type == "pass") | is.null(na.action.type)) {
            	     if(length(which.na.pass)>0)  y <- y[, - which.na.pass]
           	}

        	if (!is.null(w)) {
            	sum.y <- t(y[wind,,drop=FALSE]) %*% matrix(1,ncol=1,nrow=n)
        	} else sum.y <- t(y) %*% matrix(1,ncol=1,nrow=n)
        
		# Find abundance ranks OF MVABUND.OBJECT.1.
        	var.subset <- order(sum.y, decreasing = TRUE)
        
        	# Ensure no more than n.vars in var.subset.
        	if (n.vars<length(var.subset)) var.subset<-var.subset[1:n.vars]               	
	# Arrange data to plot requested var.subset (default - n.vars most abund).
        } else if (p<max(var.subset)) {
          	stop ("You have passed an invalid var.subset")
        # Do some dimension checks for the subset.
        } else if (n.vars!=var.subset.dim) { 
		n.vars<- var.subset.dim
        }   

    	r       <- r[,var.subset, drop=FALSE]
    	r.orig  <- r.orig[,var.subset, drop=FALSE]
    	yh      <- yh[,var.subset, drop=FALSE]
    	n.vars  <- ncol(r)

    	########## END edit var.subset, n.vars and r & fitted values  ##############

    	var.names <- colnames(r)
    	if(is.null(var.names)) var.names <- as.character(1:n.vars)

      if( "col" %in% names(dots) )
        col <- dots$col
      else
        col = rainbow(n.vars+1)[2:(n.vars+1)]

    	#################### BEGIN get window dimensions  ##########################
  
    	if(!is.null(mfcol)) {
		mfrow <- mfcol  
	}
    
        # Get all the graphical parameters.
        opp <- par("col.main","mfrow","mfcol","oma")

        if (length(mfrow)==1){
        	# i.e. mfrow is an integer either the default or a passed value,
        	# ie calc nrows & ncols
        	if ((overlay & write.plot=="show" & mfrow <5) | (mfrow <4)) {
            	    if(write.plot=="show" & is.null(dev)) {
               		if (mfrow==1) { 
				height <- 14
                 		width <- 12
               		} else {
                 		width   <- 10	#MODDED BY SW
				height <- 8
               		}
                	dev.off()
                	do.call(dev.name, args=list(height=height,width=width))
                     }
            
                     par(mfrow=c(1, mfrow)) 
                     row <- 1
            	     columns <- mfrow            
		} else {
                     columns <- ceiling(sqrt(mfrow))
                     row <- columns-1
                     if (row*columns<mfrow) row <- columns            
              	     if (write.plot=="show" & is.null(dev)) { 
            		 if (columns > row){
                  	    width   <- 9.2
                            height  <- max(row*width/columns * 1.2,5)
                	 } else {
                    	    height  <- 11
                    	    width   <- max(height*columns/row * 0.83,4)
                	 }
                	
			dev.off()
                	do.call(dev.name, args=list(height=height,width=width))
            	     }
            	     par(mfrow=c(row, columns))
                }

        	pw <- row* columns - mfrow 
                nonwindow <- FALSE        
    	} else { 
    		if(!is.null(c(mfrow, mfcol))){
        		row <- mfrow[1]
       	                columns <- mfrow[2]
        		nonwindow <- FALSE
    		} else {
        		nonwindow <- TRUE  
        		row     <- opp$mfrow[1]
        		columns <- opp$mfrow[2]
    		}
    
    		if(write.plot=="show" & is.null(dev)) {  
        		if (columns > row){
            		    width   <- 16
            	            height  <- max(row*width/columns*1.2,5)
        		} else {
            		    height  <- 11
            		    width   <- max(height*columns/row*0.83,4)
        		}
        
			#MODDED by SW - Add feature for single plot for non-overlay
			if (length(which)==1){
			    width <-8
			    height <-10
		            mfrow <- c(4,3)
			}
			dev.off()
        		do.call(dev.name, args=list(height=height,width=width))
    		}
    

    		if (any(mfrow!=par("mfrow"))) {
    			par(mfrow=mfrow) 
		}
    	
		if (!is.null(c(mfrow, mfcol))) {
			mfrow <- row* columns
		}
    		
		pw <- 0
    	} 
 
    	if (!is.null(mfcol)) { 
        	par(mfcol=c(row,columns))   
    	} else if (!is.null(mfrow)) {
        	par(mfrow=c(row,columns))   
	}

	if (length(which)==1){
		t <- ceiling(min(n.vars,12)/3)
		par(mfrow=c(t,3))
	}
		

    	one.fig <- prod(par("mfrow")) == 1  # ie if mfrow=NULL
    

    	if (is.null(sub.caption) ) {        
		# construct the sub.caption
		fktName <- "manylm"
        	
		terms <- deparse(x$terms, width.cutoff = 70)[1] 
        	nc <- nchar(terms)
        	if (length(x$terms)>1 | nc>60) 
			terms <- paste(substr(terms,1,min(60, nc)), "...") 
        	sub.caption <- paste(fktName, "(", terms,")", sep="")
    	}

    	if(!is.null(sub.caption) && !one.fig) { 
        	oma <- par("oma")
      	if (oma[3] <2 & (is.null(dev) |  !is.null(mfrow) | !is.null(mfcol))) {
            	oma[3]<- 5
            	par(oma = oma)
        	}   
        }
    	dr <- par("oma")[3]!=0
    	
	# Ensure that mfrow = NULL for the last command.
    	if (is.null(mfrow)) {mfrow <- row* columns}
    	if(all(opp$mfrow == c(row,columns))) opp$mfrow <- opp$mfcol <- NULL
    	if(keep.window & write.plot=="show") opp$mfrow <- opp$mfcol <- opp$oma <- NULL
    	
	# Upon exiting the function, reset all graphical parameters to its value
    	# at the beginning.
    	if(write.plot=="show") on.exit( par(opp), add=TRUE )

    	########################### END get window dimensions  #####################


    	########################### BEGIN selection of colors ######################

    	lcols <- length(col)
        
    	if (lcols==p & lcols != n.vars) { 
		# Adjust the colors to the subset
        	col <- col[var.subset]
        	lcols <- n.vars 
    	} else if (lcols>n.vars) {
    		col  <- col[var.subset]
    		lcols <- n.vars
    		warning("Only the first ", n.vars, " colors will be used for plotting.")
    	} else if (lcols>1 & lcols<n.vars) {
		col <- col[1]
    	    	lcols <- 1
    		warning("The vector of colors has inappropriate length. 
				Only the first color will be used")
    	}
    
	color <- col    
    	if (lcols == 1) color <- rep(color, times = n.vars) 
    	###########################  END selection of colors  ######################
    	if (any(show[2:4])) { 
            if (df.residual(x)==0 && studentized && any(show[2:3])) 
            	stop("Plot(s) ", c(2,3)[show[2:3]], " cannot be drawn: standardized residuals cannot be calculated, as there are no degrees of freedom")
            else {  
               	wr <- na.omit( as.matrix(weighted.residuals(x)) )[,var.subset,drop=FALSE]                # variance <- (t(wr) %*% wr) / df.residual(x)
            } 
            ######## to check: use only the diagonal of this variance 
            ######## (as in deviance.mlm,  rep(1, nrow(wr)) %*% wr^2)
       	    s <- sqrt(deviance(x)/df.residual(x))  
            s <- s[var.subset]
#            if(any(na.action.type == "pass") | is.null(na.action.type)) {
#                  s <- s[- which.na.pass] }
            if (show[4]) {          
                # the non-weighed residuals are used!
                cook <- cooks.distance.manylm(x, sd = s, res = r.orig)
                # for lm: the cooks distance is calculated for var.subset
               	# because of s and r from varsubset
               	# for glm: the cooks distance is calculated for all variables,
               	# select var.subset afterwards
                if(any(na.action.type == "exclude"))  cook <- cook[ - na.action , ]
            }
            if (any(show[c(2:3)])) {
	       	r.w <- if (is.null(w)) { r } # weighed residuals
                       else sqrt(w) * r
            	if(studentized) {
		     ylab23 <- "Standardized residuals"
                     if(any(na.action.type == "exclude"))
                        hii <- diag(x$hat.X[-na.action])
                     else
                        hii <- diag(x$hat.X)
                     stud  <- (1- hii)^(-(1/2))
                     rs <- (diag(stud) %*% r.w)/s  
                } else {
	             ylab23 <- "Residuals"
                     rs <- r.w
                }
       	        rs[is.infinite(rs)] <- NaN
           }
        }
    
	if (any(show[c(1, 3)])) l.fit <- "Fitted values"
    	if (is.null(id.n)) id.n <- 0
    	else {
      	    id.n <- as.integer(id.n)
            if (id.n < 0 || id.n > n) 
            	stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)

    	    if (id.n > 0) {
        	if (is.null(labels.id)) labels.id <- paste(1:n)
        	iid <- 1:id.n        
           # Obtain vector of positions of the abs. highest values, use with a vector.
        	if (overlay) {
                   show.r <- matrix(ncol=n.vars, nrow=id.n) 
                   for (i in 1:n.vars) {
              	      show.r[,i] <- (order(abs(r[,i]), decreasing = TRUE))[iid] +(i-1)*n
            	   }
            	   show.r <- c(show.r) 
        	} else {
            	   show.r <- matrix(ncol=n.vars, nrow=n) 
            	   for (i in 1:n.vars) {
			show.r[,i] <- (order(abs(r[,i]), decreasing = TRUE))
		   }
            	   show.r <- show.r[iid,]
        	}

        	if (any(show[2:3])) {
            	   if (overlay) {
                	show.rs <- matrix(ncol=n.vars, nrow=id.n)
                	for (i in 1:n.vars) {
                	show.rs[,i]<-(order(abs(rs[,i]),decreasing = TRUE))[iid] +(i-1)*n
                	}
                	show.rs <- c(show.rs)
                   } else {
                        show.rs <- matrix(ncol=n.vars, nrow=n) 
                        for (i in 1:n.vars) {
			    show.rs[,i] <- (order(abs(rs[,i]), decreasing = TRUE))
                        }
                        show.rs <- show.rs[iid,]
            	   }
        	}
       		
		##### what on earth is THIS? #####
 		text.id <- function(x, y, labels, adj.x = TRUE, col="black") {
        	    # function to write a text at a plot at the position labpos
                    	labpos <- if (adj.x) label.pos[1 + as.numeric(x > mean(range(x)))] else 3
                        text(x, y, labels, cex = cex.id, xpd = TRUE,pos = labpos, offset = 0.25, col=col)
        	}
    	   }
       }

	#######THIS IS SECTION IS FOR OVERLAY=TRUE #########
	# this creates only one plot per diagnostics	   #
	####################################################
 
    	if (overlay | n.vars==1) {      
		# plot all variables together
        	if (missing(ask)) {
			ask <- ( dev.interactive() &
          				((prod(mfrow) < length(which)) | 
						(nonwindow & !is.null(dev)) ) )
		}

        	if (ask) {   
            	op <- par(ask = TRUE)       # if TRUE, this should be preserved
           		on.exit(par(op), add=TRUE )     
		} 
        
		if(substr(legend.pos, 1,1)!="none") {   # add a legend
    
            	ncoll<- ceiling(n.vars/(50/(row+0.5*row^2))) 
            	cexl<-0.9
            	if (ncoll>3) {
				ncoll<-3	
				cexl <- 0.6
			}
            	leg <- substr(var.names, 1,(8/ncoll)+1)
        	}               

		#SW - Reset mfrow to be approprite for one plot, set legend position
		if (length(which) > 1) {
	
			if (length(which)==2) {
				par(mfrow=c(1,2),oma=c(0,1,2,5))

			} else if (length(which)==3) {
				par(mfrow=c(2,2),oma=c(0,1,2,5))

			} else {
				par(mfrow=c(2,2),oma=c(0,1,2,5))
			}

		} else {
			par(mfrow = c(1,1), oma=c(0,1,2,5))
			legend.pos="right"
		}
		
		if (show[1]) {
                   xlim <- range(yh)
                   ylim <- range(r, na.rm = TRUE)
                   if (id.n > 0) 
                      # for compatibility with R 2.2.1
                      ylim <- ylim + c(-0.08, 0.08) * diff(ylim)
                      do.call( "plot", c(list( t(yh), t(r), xlab = l.fit, ylab = "Residuals", main = main, ylim = ylim, xlim=xlim, type = "n"), dots))
          
		# Use vector built of transposed x bzw y in order to plot
            	# in the right colors.
            	      do.call( "panel", c(list(t(yh), t(r), col=color), dots))
            
		      if (one.fig) 
                      do.call( "title", c(list(sub = sub.caption), dots))
                      mtext(caption[1], 3, 0.25, col=colmain, cex=cex.caption) # the title
            
                      if (id.n > 0) {
            		# add id.n labels
                    	y.id <- (c(r))[show.r]
                    	y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
                    	text.id( (c(yh))[show.r], y.id, (rep(labels.id,              				times=n.vars))[show.r], col=rep(col, each=id.n))
                      }

                      abline(h = 0, lty = 3, col = "grey")
            	
		      if(substr(legend.pos, 1,1)[1]!="n"){    
			# add a legend
                	legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.15,xpd=NA)
                      }
        	}

        
		if (show[2]) {
            	ylim <- range(rs, na.rm = TRUE)
            	ylim[2] <- ylim[2] + diff(ylim) * 0.075
            
            	qq <- do.call( "qqnorm", c(list(t(rs), main = main, ylab = ylab23,
						col=color, asp=1), dots))
            
			if (qqline) qqline(t(rs), lty = 3, col = "gray50")
            
			# Use vector built of transposed x bzw y in order to plot
            	# in the right colors.
            
			if (one.fig) do.call( "title", c(list(sub = sub.caption), dots))
            
			mtext(caption[2], 3, 0.25, col=colmain, cex=cex.caption) # the title
            
			if (id.n > 0)                    # add id.n labels  
                		text.id(qq$x[show.rs], qq$y[show.rs], (rep(labels.id,
                  				each=n.vars))[show.rs], col=rep(col, each=id.n))
            
			if(substr(legend.pos, 1,1)[1]!="n"){   
				# add a legend
                		legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.15,xpd=NA)  
			}
        	}

        	if (show[3]) {
            	sqrtabsr <- sqrt(abs(rs))
            	ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
            	yl <- as.expression(substitute(sqrt(abs(YL)),list(YL = as.name(ylab23))))
            
		yhn0 <- yh

            	do.call( "plot", c(list(t(yhn0), t(sqrtabsr), xlab = l.fit,
              				ylab = yl, main = main, ylim = ylim, type = "n"), dots))
            
			# Use vector built of transposed x bzw y in order to plot
            	# in the right colors.
            	do.call( "panel", c(list(t(yhn0), t(sqrtabsr), col=color), dots))
            
			if (one.fig) 
                		do.call( "title", c(list(sub = sub.caption), dots))
            	mtext(caption[3], 3, 0.25, col=colmain, cex=cex.caption)
            
			if (id.n > 0) 
                		text.id(yhn0[show.rs], sqrtabsr[show.rs], (rep(labels.id,
                  				each=n.vars))[show.rs], col=rep(col, each=id.n))
            
			if(substr(legend.pos, 1,1)[1]!="n"){   
				# add a legend
                		legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.15,xpd=NA)
             	}
        	}

        	if (show[4]) {
            	ymx <- max(cook, na.rm = TRUE)
            	if (id.n > 0) {
                		show.r <- matrix(ncol=n.vars, nrow=id.n) 
                
				for (i in 1:n.vars) {
                    		show.r[,i] <- (order(-cook[,i]))[iid] +(i-1)*n
                		}
                
				show.r <- c(show.r) 
                		ymx <- ymx * 1.075
            	}
            
			obsno <- rep(1:n, each = n.vars) + rep.int((1:n.vars)/(2*n.vars),times =n)
               
            	# Use vector built of transposed x bzw y in order to plot
            	# in the right colors.
            	do.call( "plot", c(list(obsno, c(t(cook)), xlab = "Obs. number",
                			ylab = "Cook's distance", main = main, ylim = c(0, ymx), 
						type = "h",  col=color), dots))      
            
			if (one.fig) 
                		do.call( "title", c(list(sub = sub.caption), dots))
            	mtext(caption[4], 3, 0.25, col=colmain, cex=cex.caption)
            
			if (id.n > 0) {
                		txtxshow <- show.r + rep((1:n.vars)/(2*n.vars), each =id.n)-
                 			 		rep((0:(n.vars-1))*n,  each =id.n)
                
                		text.id(txtxshow, (c(cook))[show.r], (rep(labels.id,times=n.vars))[show.r], 
							adj.x = FALSE, col=rep(col, each=id.n))
            	}

            	if(substr(legend.pos, 1,1)[1]!="n"){   
				# add a legend
                		legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.15,xpd=NA)  
			}
		}

        	if (legend.pos=="nextplot" ){   
			# add a legend
            	#plot(0,type="n", xlab="",ylab="", axes=FALSE)
            	#legend("left", legend=leg, col=color,pch=1, ncol=ncoll, cex=cexl) 		#MODDED BY SW 
			legend("right", legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.35,xpd=NA)
		}
            
        	# add a subcaption
        	if (!one.fig && !is.null(sub.caption) && dr) {  
            	mtext(sub.caption, outer = TRUE, cex = 1.1*par("cex.main"),col= par("col.main") )
        	}
        
       	if(n.vars < p) {
      		if(miss.varsubset) 
                   tmp <- " \n(the variables with highest total abundance)"   
		else 
            	   tmp <- " (user selected)"  
		message("Only the variables ",paste(colnames(r), collapse = ", "), " were included in the plot", tmp, ".")
      	}

        	return(invisible())

      #######THIS IS SECTION IS FOR OVERLAY=FALSE #########
	# creates a set of diagnostic plots for each spp    #
	#####################################################  
    	} else {		
    
    		nplots <- length(which)*n.vars 
    
		if (missing(ask)) 
			ask <- ( dev.interactive() & ((mfrow < nplots )|(nonwindow & !is.null(dev)) ) )
    
		if (ask) {    
			op <- par(ask = TRUE)
           		on.exit(par(op), add=TRUE )     
		}
    
    		if (!one.fig && dr) {
    			# Define a function 'scaption' to draw the sub.caption and / or open a new
    			# window after mfrow plots.
    			# par("oma")[3]: the size of the outer margins of the top in lines of text.
            	scaption <- function(i) { 
                		if (i==mfrow) {
                			mtext( sub.caption, outer = TRUE, cex =1.1*par("cex.main"),col= par("col.main"))
                			k <- 0
                
					while(k<pw) {
			                k<-k+1
                			}
                			return(1) 
				} else return(i+1)
        		}

    		} else scaption<- function(i) {} 
        
    		scapt <- 1
	
		for (i in 1:n.vars){  
                    xlim <- range(yh)
		# draw plots for all variables
       		    if (show[1]) {
            		ri <- r[,i]
            		yhi <- yh[,i]
            		ylim <- range(ri, na.rm = TRUE)
            
			if (id.n > 0) 
                	# for compatibility with R 2.2.1
                           ylim <- ylim + c(-0.08, 0.08) * diff(ylim)
                
            		do.call( "plot", c(list(yhi, ri, xlab = l.fit, ylab = "Residuals", main = main, ylim = ylim, xlim=xlim, type = "n", asp=asp), dots))
            
            		do.call( "panel", c(list(yhi, ri, col=color[i]), dots))
            
			if (one.fig) 
             		do.call( "title", c(list(sub = sub.caption), dots))
            		if (missing(caption)) 
                            capt <- paste(var.names[i], caption[1], sep="\n") 
             		else capt <- caption[1]
            
			mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
			if (id.n > 0) {                
			# draw id.n labels in the plot
                	   y.id <- ri[show.r[,i]]
                	   y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
       			   text.id(yhi[show.r[,i]], y.id, labels.id[show.r[,i]])
            		}
            
			abline(h = 0, lty = 3, col = "grey")
            		scapt <- scaption(scapt)
          	}
        		
		if (show[2]) {
            		rsi <- rs[,i]
            		ylim <- range(rsi, na.rm = TRUE)
            		ylim[2] <- ylim[2] + diff(ylim) * 0.075
            
            		qq <- do.call( "qqnorm", c(list(rsi, main = main, ylab = ylab23,
              						ylim = ylim, col=color[i], asp=asp), dots))
                
            		if (qqline) 
                		qqline(rsi, lty = 3, col = "gray50")
            
				if (one.fig) 
                			do.call( "title", c(list(sub = sub.caption), dots))
            
				if (missing(caption)) capt <- paste(var.names[i], caption[2], sep="\n") 
                		else capt <- caption[2]
            
				mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
				if (id.n > 0) text.id(qq$x[show.rs[,i]], qq$y[show.rs[,i]],labels.id[show.rs[,i]])
            		scapt <- scaption(scapt)
        		}
	
        		if (show[3]) {
            		sqrtabsr <- sqrt(abs(rs[,i]))
            		ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
            		yl <- as.expression( substitute(sqrt(abs(YL)),list(YL = as.name(ylab23))))
            
            		yhn0 <- yh[,i]

            		do.call( "plot", c(list(yhn0, sqrtabsr, xlab = l.fit, ylab = yl,
                				main = main, ylim = ylim, type = "n"), dots)) 
            
				do.call( "panel", c(list(yhn0, sqrtabsr, col=color[i]), dots) )
            
				if (one.fig) 
                			do.call( "title", c(list(sub = sub.caption), dots))
            
				if (missing(caption)) capt <- paste(var.names[i], caption[3], sep="\n") 
                 		else capt <- caption[3]
            
				mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
				if (id.n > 0)                 # draw id.n labels in the plot
                			text.id(yhn0[show.rs[,i]], sqrtabsr[show.rs[,i]],labels.id[show.rs[,i]] )
            
				scapt <- scaption(scapt)
        		}
        
			if (show[4]) {
            		if (id.n > 0) {
                			show.r4 <- order(-cook[,i])[iid]
                			ymx <- cook[show.r4[1],i] * 1.075
            		} else ymx <- max(cook[,i], na.rm = TRUE)

            		do.call( "plot", c(list(cook[,i], xlab = "Obs. number",
                				ylab = "Cook's distance", main = main, ylim = c(0, ymx), 
							type = "h", col=color[i]), dots))    
            		if (one.fig) 
                			do.call( "title", c(list(sub = sub.caption), dots))
            
				if (missing(caption)) capt <- paste(var.names[i], caption[4], sep="\n") 
            		else capt <- caption[4]    
            
				mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
				if (id.n > 0) {               # draw id.n labels in the plot
                			text.id(show.r4, cook[show.r4,i], labels.id[show.r4], adj.x = FALSE)
            		}
            
				scapt <- scaption(scapt)    
        		}
     		}
        
     		if(n.vars < p) {
     			if(miss.varsubset) tmp <- " \n(the variables with highest total abundance)"   
			else {
          			tmp <- " (user selected)" 
			}
         
			message("Only the variables ", paste(colnames(r), collapse = ", "),
        					" were included in the plot", tmp, ".")
      	}

      	return(invisible())
	}
}
 
