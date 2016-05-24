################################################################################
# PLOT.MANYLM, PLOT.MANYGLM:                                                   #
# Plot for evaluation of goodness of fit for lm.mvabund objects                #
################################################################################

default.plot.manyglm  <- function(x, which = 1, res.type="pit.norm", caption = c("Residuals vs Fitted","Normal Q-Q", "Scale-Location", "Cook's distance"), overlay=TRUE, n.vars=Inf, var.subset=NULL, panel = if (add.smooth) panel.smooth else points, sub.caption = NULL, main = "", ask, ..., id.n = if (overlay) 0 else 3, labels.id=rownames(x$Pearson.residuals), cex.id = 0.75, qqline = TRUE, add.smooth = if(!is.null(getOption("add.smooth"))){ getOption("add.smooth") } else TRUE, label.pos = c(4, 2), cex.caption=1.5, asp = 1, legend.pos= "nextplot",	mfrow= if(overlay) {length(which)+(legend.pos=="nextplot")} else if(write.plot=="show") c(min(n.vars,3),length(which)) else length(which), mfcol=NULL, write.plot="show", filename="plot.mvabund", keep.window= if(is.null(c(mfrow,mfcol))) TRUE else FALSE, legend=FALSE) 
{        
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if ("cex" %in% names(dots)) cex <- dots$cex
    else cex <- 1.5
    if ("cex.lab" %in% names(dots)) clab <- dots$cex.lab
    else clab <- 1.5
    if ("col.lab" %in% names(dots)) colab <- dots$col.lab
    else colab <- par("col.lab")
    if ("lwd" %in% names(dots)) lwd <- dots$lwd
    else lwd  <- 2
    if ("cex.axis" %in% names(dots)) caxis <- dots$cex.axis
    else caxis <- 1.5

    dev <- dev.list()
    dev.name <- getOption("device")
    
    if (write.plot!="show") {
    	if (write.plot=="eps" | write.plot=="postscript") 
           postscript(paste(filename,".eps", sep="") )
        else if (write.plot=="pdf") 
	   pdf(paste(filename,".pdf", sep="") )
        else if (write.plot=="jpeg" )
	   jpeg(paste(filename,".jpeg", sep=""))
        else if (write.plot=="bmp" )
	   bmp(paste(filename,".bmp", sep=""))
        else if (write.plot=="png" )
	   png(paste(filename,".png", sep=""))
    	on.exit( dev.off() )
   }

   if (length(dots)>0)
   {  
    	# in the plot function.
    	deactive <- c("xlab", "ylab", "ylim", "sub", "type") 
    	deactivate <- (1:length(dots))[names(dots) %in% deactive ]
  
	    for (i in length(deactivate):1) 
   	    dots[ deactivate[i] ]<-NULL #fixed up [[]], due to compile error (v2.10).

    	dots <- lapply( dots, eval, parent.frame() )
	    if( "col.main" %in% names(dots) ) colmain <- dots$col.main
	    else colmain <- par("col.main")
   }
   else colmain <- par("col.main")

   if (!inherits(x, c("manylm", "manyglm"))) warning("use 'plot.manylm' only with \"manylm\" or \"manyglm\" objects")
        
   if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
    	stop("'which' must be in 1:4")     
    
   isGlm <- inherits(x, "manyglm")
   show <- rep.int(FALSE, times = 4)
   show[which] <- TRUE
    	
   if(ncol(x$x) > 0 ) empty <- FALSE 
   else empty <- TRUE

#   if(empty && show[4]) {
   if(show[4])
   {
     	if (length(which)==1) stop("Plot no. 4 cannot be drawn, as Cooks distance cannot be calculated for an empty model")
      else
      {
        warning("Plot no. 4 cannot be drawn, as Cooks distance cannot be calculated for an empty model")
	      show[4] <- FALSE
	      which   <- which[-which[which==4]]
    	}
   }

   if (substr(res.type,1,3)=="pit")
   {
        r <- as.matrix(x$PIT.residuals) 
        if (res.type=="pit.norm") r <- residuals(x)
   }
   else r <- as.matrix(x$Pearson.residuals) # residuals(x)
   
   yh <- as.matrix(x$linear.predictor)

#    yh <- as.matrix(x$fitted.values)
#    w <- x$sqrt.weight * x$sqrt.weight
   w <- NULL

    # Change logical var.subset to numerical var.subset, if necessary. Note that NA values are logical as well, but should be excluded here.
    if(!is.null(var.subset) & !is.numeric(var.subset))
      	var.subset <- which(var.subset[!is.na(var.subset)])        
#    miss.varsubset<-!is.numeric(var.subset) # If this function is called within another, the missing function could be tricked out.

    if(is.null(labels.id)) labels.id <- as.character(1:nrow(r))

    if (!is.null(w)) {
       	wind <- w != 0
	if (isGlm & is.matrix(w)){
            wind  <- rowSums( wind) != 0
            w     <- w[wind,, drop=FALSE]
        } else w <- w[wind]        	

       	r <- r[wind,, drop=FALSE]
       	yh <- yh[wind,, drop=FALSE]  
       	labels.id <- labels.id[wind]
    }
    n <- nrow(r)
    p <- ncol(r)

    ######## BEGIN edit var.subset, n.vars and r & fitted values  #########
    # subset allows double variables
    # Do some dimension checks for the subset.
    if (missing(var.subset) | is.null(var.subset) | !is.numeric(var.subset)) {
       # Plot the n.var variables with highest abundances
       if ( p < n.vars ) {
   #       warning("You have passed an invalid number of variables 'n.vars' to be included in the plot. All variables will be included instead.")        
          n.vars <- p
       }      
       y <- as.matrix(x$y)
       if (!is.null(w)) 
          sum.y <- t(y[wind,,drop=FALSE]) %*% matrix(1,ncol=1,nrow=n)
       else sum.y <- t(y) %*% matrix(1,ncol=1,nrow=n)
        
       # Find abundance ranks OF MVABUND.OBJECT.1.
       var.subset <- order(sum.y, decreasing = TRUE)
       typeofvarsubset <- " \n(the variables with highest total abundance)"   
    }       
    else {  # if var.subset is specified
       if ( p < max(var.subset) ) 
          stop ("You have passed an invalid var.subset")
       var.subset.dim <- length(var.subset) 

       if ( missing(n.vars) | n.vars != var.subset.dim ) {
          n.vars <- var.subset.dim
	  warning("Number of variables 'n.var' is set to the length of 'var.subset'.")
       }
       typeofvarsubset <- " (user selected)"   
   } 
   
   ############# Extract relevant data ###################
   r <- r[,var.subset, drop=FALSE]
   yh <- yh[,var.subset, drop=FALSE]
   w <- w[,var.subset, drop=FALSE]
   ######### END edit var.subset, n.vars and r & fitted values  ###########

   var.names <- colnames(r)
   if(is.null(var.names)) var.names <- as.character(1:n.vars)
    
   ### SET COLORS AND GET SOME GRAPHICS PARAMETERS
   # Upon exiting the function, reset all graphical parameters to its value
   # at the beginning.
   if(!is.null(mfcol)) mfrow <- mfcol      
   # Get all the graphical parameters.
   opp <- par("col.main","mfrow","mfcol","oma")
   if( "col" %in% names(dots) )
   {
     col <- dots$col[var.subset]
     dots$col=NULL
   } 
   else
     col = rainbow(n.vars+1)[2:(n.vars+1)]
   if (write.plot=="show")
     on.exit( par(opp), add=TRUE ) 

   ################# BEGIN get window dimensions  #########################

   if (length(mfrow)==1){
      # i.e. mfrow is an integer either the default or a passed value,
      # ie calc nrows & ncols
      if ((overlay & write.plot=="show" & mfrow <5) | (mfrow <4)) {
          if(write.plot=="show" & is.null(dev)) {
              if (mfrow==1) { 
	           height <- 14
                   width <- 12
              } else {
                   width   <- 10 #MODDED BY SW
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
            if(write.plot=="show" & is.null(dev)) { 
            	if (columns > row){
                   width <- 9.2
                   height <- max(row*width/columns * 1.2,5)
       		} else {
		   height <- 11
		   width  <- max(height*columns/row * 0.83,4)
		}
                dev.off()
		do.call(dev.name, args=list(height=height,width=width))
            }
            par(mfrow=c(row, columns))
        }
       	pw <- row* columns - mfrow 
        nonwindow <- FALSE        
    } else { # if length(mfrow)==1)
        if(!is.null(c(mfrow, mfcol))){
	   row <- mfrow[1]
	   columns <- mfrow[2]
	   nonwindow <- FALSE
	} else {
	   nonwindow <- TRUE  
	   row  <- opp$mfrow[1]
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
              height <- 10
              mfrow <- c(1,1)
           }
           dev.off()
           do.call(dev.name, args=list(height=height,width=width))
        }
        if (any(mfrow!=par("mfrow"))) par(mfrow=mfrow) 		   	
        if (!is.null(c(mfrow, mfcol))) mfrow <- row* columns		    
	pw <- 0
    } 
 
    if (!is.null(mfcol)) par(mfcol=c(row,columns))   
    else if (!is.null(mfrow)) par(mfrow=c(row,columns)) 
    if (length(which)==1){
 	t <- ceiling(min(n.vars,12)/3)
#	par(mfrow=c(t,3))
        par(mfrow=c(1,1))
    }		
    
    one.fig <- prod(par("mfrow")) == 1  # ie if mfrow=NULL
    
    if (is.null(sub.caption) ) {   
        # construct the sub.caption
       	fktName <- "manyglm" 
	terms <- deparse(x$terms, width.cutoff = 70)[1] 
       	nc <- nchar(terms)
       	if (length(x$terms)>1 | nc>60) 
           terms <- paste(substr(terms,1,min(60, nc)), "...") 
        sub.caption <- paste(fktName, "(", terms,")", sep="")
    }

    if(!is.null(sub.caption) && !one.fig) {
       	oma <- par("oma")
      	if (oma[3] < 2 & (is.null(dev) | !is.null(mfrow) | !is.null(mfcol))) {
       	    oma[3]<- 5
            par(oma = oma)
        }   
    }
    dr <- par("oma")[3]!=0
    	
    # Ensure that mfrow = NULL for the last command.
    if (is.null(mfrow)) {mfrow <- row* columns}
    if(all(opp$mfrow == c(row,columns))) opp$mfrow <- opp$mfcol <- NULL
    if(keep.window & write.plot=="show") opp$mfrow <- opp$mfcol <- opp$oma <- NULL
    	
    ##################### END get window dimensions  ##################

    ####################### BEGIN selection of colors ##################
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

    #######################  END selection of colors  #################


   if (any(show[2:3])) {
      if (df.residual(x)==0) 
       	stop("Plot(s) ", c(2,3)[show[2:3]], " cannot be drawn: standardized residuals cannot be calculated, as there are no degrees of freedom")
      else {
         #############################
	 ## Q: why weighed residuals?
	 ############################

         rs <- if (is.null(w)) { r } # weighed residuals
	       else sqrt(w) * r   
         if (substr(res.type,1,3)=="pit") {
             if (res.type=="pit.norm") ylab23<-"Dunn-Smyth Residuals"
             else ylab23 <- "PIT Residuals."
         }
         else ylab23 <- "Standard Pearson residuals."
         rs[is.infinite(rs)] <- NaN
      }
   }   

   if (any(show[c(1, 3)])) 
       # l.fit <- "Fitted values"
       l.fit <- "Linear predictor value"
        
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
              for (i in 1:n.vars) show.r[,i] <- (order(abs(r[,i]),decreasing=TRUE))[iid] +(i-1)*n
              show.r <- c(show.r) 
           } else {
	      show.r <- matrix(ncol=n.vars, nrow=n) 	      
	      for (i in 1:n.vars) show.r[,i] <- (order(abs(r[,i]), decreasing = TRUE))
	      show.r <- show.r[iid,]
	   }

           if (any(show[2:3])) {
	      if (overlay) {
	         show.rs <- matrix(ncol=n.vars, nrow=id.n)
		 for (i in 1:n.vars) show.rs[,i] <-(order(abs(rs[,i]), decreasing = TRUE))[iid] +(i-1)*n
		 show.rs <- c(show.rs)
              } else {
	         show.rs <- matrix(ncol=n.vars, nrow=n) 
		 for (i in 1:n.vars) show.rs[,i] <- (order(abs(rs[,i]), decreasing = TRUE))
		 show.rs <- show.rs[iid,]
              }
          }
	  
	  ##### what on earth is THIS? #####
 	  text.id <- function(x, y, labels, adj.x = TRUE, col="black") {
          # function to write a text at a plot at the position labpos
         	labpos <- if (adj.x) 
		              label.pos[1 + as.numeric(x > mean(range(x)))]           		    	 else 3
            	text(x, y, labels, cex = cex.id, xpd = TRUE,pos = labpos, offset = 0.25, col=col)
          }
       }
   }

	#######THIS IS SECTION IS FOR OVERLAY=TRUE #########
	# this creates only one plot per diagnostics	   #
	####################################################
 
   if (overlay | n.vars==1) {      
    # plot all variables together
       if (missing(ask)) 
           ask <- ( dev.interactive() &	((prod(mfrow) < length(which)) | (nonwindow & !is.null(dev)) ) )

       if (ask) {   
          op <- par(ask = TRUE)       # if TRUE, this should be preserved
	  on.exit(par(op), add=TRUE )     
       }        
       if (substr(legend.pos, 1,1)!="none") {   # add a legend
          ncoll<- ceiling(n.vars/(50/(row+0.5*row^2))) 
	  cexl<- 1.5 #0.9
	  if (ncoll>3) {
	      ncoll<-3	
	      cexl <- 0.6
	  }
	  leg <- substr(var.names, 1,(8/ncoll)+1)
       }               
      #SW - Reset mfrow to be approprite for one plot, set legend position    
    if (length(which)==1) {
#          dev.off()
	  if (legend == TRUE) {
#             dev.new(height=6, width=8) # added for smaller window size    
             par(mfrow = c(1,1), oma=c(.5,.5,.5,4.5), mar=c(6, 4.5, 2, 5))
	     legend.pos="right"
	  }   
	  else { 
#             dev.new(height=6, width=6) # added for smaller window size    
             par(mfrow = c(1,1), oma=c(.5,.5,.5,.5), mar=c(6,4.5,2,.5))
	  }   
      }	  
      else if (length(which)==2) {
 #         dev.off()
#	  dev.new(height=6, width=12)
          par(mfrow=c(1,2),oma=c(0.5,0.5,1,10), mar=c(4, 4.5, 2, 2))
      }	  
      else if (length(which)==3) {
#          dev.new(height=12, widht=12)
          par(mfrow=c(2,2),oma=c(2,2,2,2), mar=c(4, 4, 3, 3))
      }	

      # The residual vs. fitted value plot	
      yhtmp <- c(yh)
      yh.is.zero <- yhtmp < (-6)
#       yh.is.zero <- yhtmp < max(-6,(-max(yhtmp)))#this line is wrong - it kicks out any value more negative than max(yh)
       yh0 <- yhtmp[!yh.is.zero]
      xlim <- range(yh0)
      if (id.n > 0) # for compatibility with R2.2.1
          ylim <- ylim + c(-0.08, 0.08) * diff(ylim) 
          
      # drop small values in the response
      if (show[1]) {
          # Use vector built of transposed x bzw y to plot in the right color
#	  yi.is.zero <- (yh[,1]<(-9)) # log(1e-4)
#	  plot(x=t(yh[!yi.is.zero,1]), y=t(r[!yi.is.zero,1]),type="p",col=palette()[1], ylab = "Pearson residuals", xlab=l.fit, main = main, ylim=ylim, xlim=xlim, cex.lab=clab, cex.axis=caxis, cex=cex, lwd=lwd, font.main=2)	
#          for (i in 2:n.vars) {	  
#              yi.is.zero <- (yh[,i] < (-9)) # log(1e-4)
#              points(x=t(yh[!yi.is.zero,i]), y=t(r[!yi.is.zero,i]),type="p",col=palette()[i], cex=cex, lwd=lwd)
#          }	  
          
          rtmp <- c(r)
          r0 <- rtmp[!yh.is.zero]
#          ylim <- range(max(range(abs(r0), finite = TRUE)) * c(-1, 1), na.rm = TRUE) #DW, 10/02/15: to make ylims symmetric about zero
          ylim <- max( range(abs(r0), finite = TRUE, na.rm=TRUE) ) * c(-1, 1) #DW, 10/02/15: to make ylims symmetric about zero
          # DW, 21/01/16: use of range suggested by Eduard Szocs to remove Inf values
          
          colortmp <- rep(color, each=n)
          color0 <- colortmp[!yh.is.zero]
            
          if (substr(res.type,1,3)=="pit") {
              if (res.type=="pit.norm") ylab="Dunn-Smyth Residuals"
              else ylab="PIT Residuals"
          }
          else ylab="Pearson residuals"
	  plot(yh0, r0, xlab = l.fit, ylab = ylab, main=main, ylim=ylim, xlim=xlim, font.main=2, col=color0, cex.lab=clab, cex.axis=caxis, cex=cex, lwd=lwd)

          # Add sub.caption, e.g, manyglm(tasm.cop ~ treatment)
          if (one.fig) 
              title(sub = sub.caption, cex.sub=cex.caption-0.1, line=4.5, font.sub=2)
	  # Add the title Residual vs Fitted 
          mtext(caption[1], 3, 0.25, col=colmain, cex=cex.caption, font=2)
	  
	  if (id.n > 0) { # add id.n labels
             y.id <- (c(r))[show.r]
             y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
             text.id( (c(yh))[show.r], y.id, (rep(labels.id, times=n.vars))[show.r], col=rep(col, each=id.n))
          }

          if (res.type=="pit.uniform") hmark=0.5 else hmark=0
          abline(h = hmark, lty = 2, col = "black", lwd=2)            	
          if(legend==TRUE & substr(legend.pos, 1,1)[1]!="n"){    
	    # add a legend
	      legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl-0.1,inset=-0.35,xpd=NA, lwd=2, lty=0, x.intersp=0.5)
	  }

#	  mtext("(a)", side = 3, cex = 2, at=-1.8, line=0.3)

      }        
      
      # The normal QQ plot
      if (show[2]) { 
#         rs.is.zero <- rs < (1e-9) #DW, 23/10/14: this seems to be an error - why would r near zero be a problem?
         rstmp <- c(rs)
#         ylim <- range(max(range(abs(rstmp), finite = TRUE)) * c(-1, 1), na.rm = TRUE) #DW, 23/10/14: to make ylims symmetric about zero
         ylim <- max( range(abs(rstmp), finite = TRUE, na.rm=TRUE) ) * c(-1, 1) #DW, 23/10/14: to make ylims symmetric about zero
         #DW, 21/1/16: finite=TRUE added as suggested by Eduard Szocs 
   	 qq <- do.call( "qqnorm", c(list(rstmp, main = main, ylab = ylab23, ylim=ylim, col=color, asp=1, cex.lab=1.5, cex=1.5, cex.axis=1.5, cex.main=1.5, lwd=2), dots))
	 if (qqline) abline(c(0,1), lty = 3, col = "gray50", lwd=2)          
         # Use vector built of transposed x bzw y in order to plot
	 # in the right colors.
	 if (one.fig) do.call( "title", c(list(sub = sub.caption), dots))
	 mtext(caption[2], 3, 0.25, col=colmain, cex=cex.caption) # the title
            
         if (id.n > 0) # add id.n labels  
	    text.id(qq$x[show.rs], qq$y[show.rs], (rep(labels.id,
	            each=n.vars))[show.rs], col=rep(col, each=id.n))

         if(legend == TRUE & substr(legend.pos, 1,1)[1]!="n"){   
	     # add a legend
	     legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl-0.1,inset=-0.35,xpd=NA, x.intersp=0.5, lwd=2, lty=0)  
	 }
      }
      # The scale vs. location plot
      if (show[3]) {
          sqrtabsr <- c(sqrt(abs(rs)))
          sqrtabsr0 <- sqrtabsr[!yh.is.zero]
	  ylim <- c(0, max(sqrtabsr0, na.rm = TRUE))
	  yl <- as.expression(substitute(sqrt(abs(YL)),list(YL = as.name(ylab23))))         

#	  yi.is.zero <- (yh[,1]<(-9))
#	  plot(t(yh[!yi.is.zero,1]), t(sqrtabsr[!yi.is.zero,1]),type="p",col=palette()[1], ylab = yl, xlab=l.fit, main = main, ylim=ylim, xlim=xlim, cex=1.5, cex.lab=1.5, cex.axis=1.5)	            
#          for (i in 2:n.vars) {	  
#              yi.is.zero <- (yh[,i] < (-9))
#              points(t(yh[!yi.is.zero,i]), t(sqrtabsr[!yi.is.zero,i]),type="p",col=palette()[i], cex=1.5, lwd=2)
#          }	  
	  plot(yh0, sqrtabsr0, xlab = l.fit, ylab=yl,
                    main = main, ylim = ylim, xlim=xlim, type = "n")
	  panel(yh0, sqrtabsr0, col=color, cex=cex, cex.lab=clab, cex.axis=caxis, lwd=lwd) 

	  if (one.fig) 
	      do.call( "title", c(list(sub = sub.caption), dots))
          mtext(caption[3], 3, 0.25, col=colmain, cex=cex.caption)
       
	  if (id.n > 0) 
	     text.id(yh[show.rs], sqrtabsr[show.rs], (rep(labels.id,
	             each=n.vars))[show.rs], col=rep(col, each=id.n))		    
#          ncoll <- ceiling(n.vars/6)
          if(legend==TRUE & substr(legend.pos, 1,1)[1]!="n") # add a legend  
	     legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl-0.1,inset=-0.35,xpd=NA, x.intersp=0.5)             
      }

      # The cook distance plot
      if (show[4]) {
#         stop("cook's distance for glm is not implemented.") 
#         ymx <- max(cook, na.rm = TRUE) # error here, what is cook?
#	 if (id.n > 0) {
#	    show.r <- matrix(ncol=n.vars, nrow=id.n) 
#	    for (i in 1:n.vars) 
#	        show.r[,i] <- (order(-cook[,i]))[iid] +(i-1)*n
#            show.r <- c(show.r) 
#            ymx <- ymx * 1.075
#         }
#         obsno <- rep(1:n, each = n.vars) + rep.int((1:n.vars)/(2*n.vars),times =n)
#         # Use vector built of transposed x bzw y in order to plot in the right colors.
#         do.call( "plot", c(list(obsno, c(t(cook)), xlab = "Obs. number", ylab = "Cook's distance", main = main, ylim = c(0, ymx), type = "h",  col=color), dots))      
#
#         if (one.fig) do.call( "title", c(list(sub = sub.caption), dots))
#         mtext(caption[4], 3, 0.25, col=colmain, cex=cex.caption)
#            
#         if (id.n > 0) {
#	    txtxshow <- show.r + rep((1:n.vars)/(2*n.vars), each =id.n)-rep((0:(n.vars-1))*n,  each =id.n)
#	    text.id(txtxshow, (c(cook))[show.r], (rep(labels.id,times=n.vars))[show.r], adj.x = FALSE, col=rep(col, each=id.n))
 #        }
 #   
 #        if(legend==TRUE & substr(legend.pos, 1,1)[1]!="n")  # add a legend
 #           legend(legend.pos, legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl,inset=-0.15,xpd=NA)
      }

      if (legend==TRUE & legend.pos=="nextplot" ) legend("right", legend=leg, col=color, pch=1, ncol=ncoll, cex=cexl-0.1,inset=-0.5,xpd=NA)
           
      # add a subcaption
      if (!one.fig && !is.null(sub.caption) && dr) 
         mtext(sub.caption, outer = TRUE, cex = 1.1*par("cex.main"),col= par("col.main") )
        
      if(n.vars < p) {
           if (missing(var.subset) | is.null(var.subset) | !is.numeric(var.subset))
	   message("Only the variables ",paste(colnames(r), collapse = ", "), " were included in the plot", typeofvarsubset, ".")
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
		while(k<pw) k<-k+1
		return(1) 
	     } else return(i+1)
          }
       } else scaption<- function(i) {} 
        
      scapt <- 1
      for (i in 1:n.vars){  
      # draw plots for all variables
         if (show[1]) {
	     ri <- r[,i]
	     yhi <- yh[,i]
	     ylim <- range(ri, na.rm = TRUE)
            
	     if (id.n > 0) 
	     # for compatibility with R 2.2.1
	         ylim <- ylim + c(-0.08, 0.08) * diff(ylim)
                
            if (res.type=="pit") ylab="Random Quantile Residuals"
            else ylab="Pearson residuals"
	    do.call( "plot", c(list(yhi, ri, xlab = l.fit, ylab=ylab, main = main, ylim = ylim, type = "n", asp=asp), dots))
            
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
        
#	    abline(h = 0, lty = 3, col = "grey")
	    scapt <- scaption(scapt)
	 }
        		
	 if (show[2]) {
	    rsi <- rs[,i]
	    ylim <- range(rsi, na.rm = TRUE)
	    ylim[2] <- ylim[2] + diff(ylim) * 0.075
	    
	    qq <- do.call( "qqnorm", c(list(rsi, main = main, ylab = ylab23, ylim = ylim, col=color[i], asp=asp), dots))

            if (qqline) 
	       qqline(rsi, lty = 3, col = "gray50")
            
	    if (one.fig) 
	       do.call( "title", c(list(sub = sub.caption), dots))
            
	    if (missing(caption)) 
	       capt <- paste(var.names[i], caption[2], sep="\n") 
	    else capt <- caption[2]

	    mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
	    if (id.n > 0) 
	       text.id(qq$x[show.rs[,i]], qq$y[show.rs[,i]],labels.id[show.rs[,i]])
            scapt <- scaption(scapt)
        }
	
        if (show[3]) {
	    sqrtabsr <- sqrt(abs(rs[,i]))
	    ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
	    yl <- as.expression( substitute(sqrt(abs(YL)),list(YL = as.name(ylab23))))
	    do.call( "plot", c(list(yh, sqrtabsr, xlab = l.fit, ylab = yl, main = main, ylim = ylim, type = "n", cex=1.5), dots)) 
            
	    do.call( "panel", c(list(yh, sqrtabsr, col=color[i]), dots) )

            if (one.fig) 
	       do.call( "title", c(list(sub = sub.caption), dots))

            if (missing(caption)) 
	       capt <- paste(var.names[i], caption[3], sep="\n") 
	    else capt <- caption[3]            

	    mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
            
#	    if (id.n > 0) # draw id.n labels in the plot
#	       text.id(yhn0[show.rs[,i]], sqrtabsr[show.rs[,i]],labels.id[show.rs[,i]] )
            scapt <- scaption(scapt)
	}
        
	if (show[4]) {
            stop("cook's distance for glm is not implemented.") 
#	    if (id.n > 0) {
#	       show.r4 <- order(-cook[,i])[iid]
#	       ymx <- cook[show.r4[1],i] * 1.075
#	    } else ymx <- max(cook[,i], na.rm = TRUE)
#
#            do.call( "plot", c(list(cook[,i], xlab = "Obs. number", ylab = "Cook's distance", main = main, ylim = c(0, ymx), type = "h", col=color[i]), dots))    
#	    if (one.fig) 
#	       do.call( "title", c(list(sub = sub.caption), dots))
#
#	    if (missing(caption)) 
#	       capt <- paste(var.names[i], caption[4], sep="\n") 
#            else capt <- caption[4]    
#            
#	    mtext(capt, 3, 0.8, col=colmain, cex=cex.caption)  # draw the title
#            
#	    if (id.n > 0)  # draw id.n labels in the plot
#	       text.id(show.r4, cook[show.r4,i], labels.id[show.r4], adj.x = FALSE)
#            scapt <- scaption(scapt)    
        }
    }  # end for

    if(n.vars < p) {
        if(missing(var.subset)|is.null(var.subset)|!is.numeric(var.subset)) 
	   tmp <- " \n(the variables with highest total abundance)"   
	else tmp <- " (user selected)" 
    }
    
    message("Only the variables ", paste(colnames(r), collapse = ", "), " were included in the plot", tmp, ".")
  }
return(invisible())

}
 
