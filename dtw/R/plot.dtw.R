###############################################################
#                                                             #
#   (c) Toni Giorgino <toni.giorgino,gmail.com>           #
#       Istituto di Ingegneria Biomedica (ISIB-CNR)                 #
#       Consiglio Nazionale delle Ricerche                           #
#       www.isib.cnr.it                                    #
#                                                             #
#   $Id: plot.dtw.R 267 2012-08-12 14:37:26Z tonig $
#                                                             #
###############################################################




## Plot a dtw non-object, switching depending on requested type

plot.dtw <- function(x, type="alignment", ...) {
  
  pt<-pmatch(type,c("alignment",
                    "twoway",
                    "threeway",
                    "density"));
  switch(pt, 	dtwPlotAlignment(x, ...),
                dtwPlotTwoWay(x, ...),
                dtwPlotThreeWay(x, ...),
                dtwPlotDensity(x, ...)
 );
}


## an alias
dtwPlot <- plot.dtw;




dtwPlotAlignment <- function(d, xlab="Query index", ylab="Reference index", plot.type="l", ...) {
  plot( d$index1,d$index2,
        xlim=c(1,d$N),ylim=c(1,d$M),
	xlab=xlab,ylab=ylab,type=plot.type,
        ...
	);
}


## Normalization plots the average cost per step instead of
## the cumulative cost

dtwPlotDensity <- function(d, normalize=FALSE,
                           xlab="Query index", ylab="Reference index", ...) {

    cm<-d$costMatrix;

    if(is.null(cm)) 
      stop("dtwPlotDensity requires dtw internals (set keep.internals=TRUE on dtw() call)");

    ## We can safely modify cm locally
    if(normalize) {
        norm <- attr(d$stepPattern,"norm");
        if(is.na(norm))
          stop("No normalization known for step pattern used");

        if(norm=="N") {
            cm <- cm / row(cm);
        } else if(norm=="N+M") {
            cm <- cm / (row(cm)+col(cm));
        } else if(norm=="M") {
            cm <- cm / col(cm);
        }
    }

    xd<-dim(cm)[1];
    yd<-dim(cm)[2];

    image(cm,col=terrain.colors(100),x=1:xd,y=1:yd,
          xlab=xlab,ylab=ylab, ...);
    contour(cm,x=1:xd,y=1:yd,add=TRUE);
    lines(d$index1,d$index2,col="blue",lwd=2);
}




## Well-known and much-copied pairwise matching

dtwPlotTwoWay <- function(d,xts=NULL,yts=NULL, offset=0,
			ts.type="l",pch=21, 
                        match.indices=NULL,
			match.col="gray70", match.lty=3,
			xlab="Index", ylab="Query value", 
			... ) {

	if(is.null(xts) || is.null(yts))  {
            xts <- d$query;
            yts <- d$reference;
        }
    
	if(is.null(xts) || is.null(yts)) 
		stop("Original timeseries are required");

        ytso<-yts+offset;

        ## pad to longest
        maxlen<-max(length(xts),length(ytso));
        length(xts)<-maxlen;
        length(ytso)<-maxlen;

	
	## save default, for resetting...
	def.par <- par(no.readonly = TRUE);

        ## make room for secondary axis, if any
        if(offset!=0) {
          par(mar=c(5,4,4,4)+.1);
        }

	## plot q+t
	matplot(cbind(xts,ytso),
                type=ts.type,pch=pch, 
                xlab=xlab, ylab=ylab,
                axes=FALSE,
                ...);

        ## box and main axis
        ## compute range covering all values
        box();
        axis(1);
        axis(2,at=pretty(xts));

        ## display secondary axis if offset
        if(offset!=0) {
          rightTicks <- pretty(yts);
          axis(4,at=rightTicks+offset,labels=rightTicks);
        }


	## plot the matching 
	# par(par.match);
        if(is.null(match.indices)) {
          ml<-length(d$index1);
          idx<-1:ml;
        } else if(length(match.indices)==1) {
          idx <- seq(from=1,
                     to=length(d$index1),
                     length.out=match.indices);
        } else {
          idx <- match.indices;
        }

	## x0, y0 	coordinates of points from which to draw.
	## x1, y1 	coordinates of points to which to draw.
	segments(d$index1[idx],xts[d$index1[idx]],
		 d$index2[idx],ytso[d$index2[idx]],
		 col=match.col,lty=match.lty);

	
	par(def.par)#- reset to default

}







## ##################################################
## Global distance density plot

# for each plot, we should set: color, width, style, type
# for match lines: color, width, style

dtwPlotThreeWay <- function(d,xts=NULL,yts=NULL,
                            type.align="l",type.ts="l",
                            match.indices=NULL,
                            margin=4, inner.margin=0.2, title.margin=1.5,
                            xlab="Query index",ylab="Reference index",main="Timeseries alignment",
                            ... ) {

     if(is.null(xts) || is.null(yts))  {
         xts <- d$query;
         yts <- d$reference;
     }

     # Sanity check
     if(is.null(xts) || is.null(yts))
       stop("Original timeseries are required");

     # Coerce to plain vectors
     xts <- as.matrix(xts);
     yts <- as.matrix(yts);

     # Verify if not multivariate
     if( ncol(xts)>1 || ncol(yts)>1 )
       stop("Only single-variate timeseries can be displayed. (You may want to extract a column for visualization purposes.)");


     def.par <- par(no.readonly = TRUE) # save default, for resetting...

     layout(matrix(c(3,1,0,2),2,2,byrow=TRUE), c(1,3), c(3,1), TRUE);

     
     imar<-inner.margin;
     
     bmar<-margin;
     lmar<-margin;
     tmar<-margin+title.margin;
     rmar<-margin;

     mlab=margin/2;
     mtex=margin/6;

     nn<-length(xts);
     mm<-length(yts);

     
     # Plot the warping function
     par(mar=c(imar,imar,tmar,rmar));

     # todo: plot over segments

     plot(d$index1,d$index2,type=type.align,
          xlim=c(1,nn),ylim=c(1,mm),
          ax=FALSE,main=main, ...
          ); # fake a diagonal, to set the axes


     # vertical match segments
     #  1 value: plot total of N elements
     if(length(match.indices)==1) {
       match.indices <- seq(from=1,
                            to=length(d$index1),
                            length.out=match.indices);
     }

     #  vector: use specified indices
     if(! is.null(match.indices) ) {      # vertical match segments
       idx <- match.indices;
       segments(d$index1[idx],0,
                d$index1[idx],d$index2[idx],
                col="grey60",lty=3);
                                        # horz.
       segments(0,d$index2[idx],
                d$index1[idx],d$index2[idx],
                col="grey60",lty=3);
     }

     
     box();

     
     # axis are 1- bot; 2- left; 3- top; 4- right
     # Plot query (horizontal, bottom)
     par(mar=c(bmar,imar,imar,rmar));

     plot(xts ~ c(1:nn), type=type.ts,
          xlab=xlab ,mgp=c(mlab,mtex,0) ,ax=FALSE,
          );
     axis(1);
     axis(2);
     box();

     # Plot reference (vertical, left)
     par(mar=c(imar,lmar,tmar,imar));

     # reverse the horiz. axis so that rotation is more natural
     plot(c(1:mm) ~ yts, xlim=rev(range(yts)), type=type.ts,
          ylab=ylab, mgp=c(mlab,mtex,0) , ax=FALSE,
          );
     axis(3);
     axis(2);
     box();

     par(def.par)#- reset to default

}

