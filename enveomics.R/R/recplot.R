enve.recplot <- structure(function(
	### Produces recruitment plots provided that BlastTab.catsbj.pl has
	### been previously executed. Requires the gplots library.
	prefix,
	### Path to the prefix of the BlastTab.catsbj.pl output files. At
	### least the files .rec and .lim must exist with this prefix.
	
	# Id. hist.
	id.min=NULL,
	### Minimum identity to be considered. By default, the minimum detected
	### identity. This value is a percentage.
	id.max=NULL,
	### Maximum identity to be considered. By default, 100.
	id.binsize=NULL,
	### Size of the identity bins (vertical histograms). By default, 0.1 for
	### identity metrics and 5 for bit score.
	id.splines=0,
	### Smoothing parameter for the splines in the identity histogram. Zero (0) for no
	### splines. A generally good value is 1/2. If non-zero, requires the stats package.
	id.metric='id',
	### Metric of identity to be used (Y-axis). It can be any unambiguous prefix
	### of "identity", "corrected identity", or "bit score".
	id.summary='sum',
	### Method used to build the identity histogram (Horizontal axis of the right panel).
	### It can be any unambiguous prefix of "sum", "average", "median", "90% lower bound",
	### "90% upper bound", "95% lower bound", and "95% upper bound". The last four options
	### correspond to the upper and lower boundaries of the 90% and 95% empirical confidence
	### intervals.
	
	# Pos. hist.
	pos.min=1,
	### Minimum (leftmost) position in the reference (concatenated) genome (in bp).
	pos.max=NULL,
	### Maximum (rightmost) position in the reference (concatenated) genome (in bp).
	### By default: Length of the genome.
	pos.binsize=1e3,
	### Size of the position bins (horizontal histograms) in bp.
	pos.splines=0,
	### Smoothing parameter for the splines in the position histogram. Zero (0) for no splines.
	### If non-zero, requires the stats package.
	
	# Rec. plot
	rec.col1='white',
	### Lightest color in the recruitment plot.
	rec.col2='black',
	### Darkest color in the recruitment plot.

	# General
	main=NULL,
	### Title of the plot.
	contig.col=grey(0.85),
	### Color of the Contig boundaries. Set to NA to ignore Contig boundaries.
	
	# Return
	ret.recplot=FALSE,
	### Indicates if the matrix of the recruitment plot is to be returned.
	ret.hist=FALSE,
	### Indicates if the vectors of the identity and position histograms are to be returned.
	ret.mode=FALSE,
	### Indicates if the mode of the identity is to be computed. It requires the modeest
	### package.
	
	# General
	id.cutoff=NULL,
	### Minimum identity to consider an alignment as "top". By default, it is 0.95 for the
	### identity metrics and 95% of the best scoring alignment for bit score.
	verbose=TRUE,
	### Indicates if the function should report the advance.
	...
	### Any additional graphic parameters to be passed to plot for all panels except the
	### recruitment plot (lower-left).
	){
   
   # Settings
   METRICS <- c('identity', 'corrected identity', 'bit score');
   SUMMARY <- c('sum', 'average', 'median', '');
   if(is.null(prefix)) stop('Parameter prefix is mandatory.');
   if(!requireNamespace("gplots", quietly=TRUE)) stop('Unavailable gplots library.');

   # Read files
   if(verbose) cat("Reading files.\n")
   rec <- read.table(paste(prefix, '.rec', sep=''), sep="\t", comment.char='', quote='');
   lim <- read.table(paste(prefix, '.lim', sep=''), sep="\t", comment.char='', quote='');

   # Configure ID summary
   id.summary <- pmatch(id.summary, SUMMARY);
   if(is.na(id.summary)) stop('Invalid identity summary.');
   if(id.summary == -1) stop('Ambiguous identity summary.');
   if(id.summary==1){
      id.summary.func <- function(x) colSums(x);
      id.summary.name <- 'sum'
   }else if(id.summary==2){
      id.summary.func <- function(x) colMeans(x);
      id.summary.name <- 'mean'
   }else if(id.summary==3){
      id.summary.func <- function(x) apply(x,2,median);
      id.summary.name <- 'median'
   }else if(id.summary==4){
      id.summary.func <- function(x) apply(x,2,quantile,probs=0.05,names=FALSE);
      id.summary.name <- '90% LB'
   }else if(id.summary==5){
      id.summary.func <- function(x) apply(x,2,quantile,probs=0.95,names=FALSE);
      id.summary.name <- '90% UB'
   }else if(id.summary==6){
      id.summary.func <- function(x) apply(x,2,quantile,probs=0.025,names=FALSE);
      id.summary.name <- '95% LB'
   }else if(id.summary==7){
      id.summary.func <- function(x) apply(x,2,quantile,probs=0.975,names=FALSE);
      id.summary.name <- '95% UB'
   }

   # Configure metrics
   id.metric <- pmatch(id.metric, METRICS);
   if(is.na(id.metric)) stop('Invalid identity metric.');
   if(id.metric == -1) stop('Ambiguous identity metric.');
   if(id.metric==1){
      id.reccol <- 3
      id.shortname <- 'Id.'
      id.fullname  <- 'Identity'
      id.units     <- '%'
      id.hallmarks <- seq(0, 100, by=5)
      if(is.null(id.max)) id.max <- 100
      if(is.null(id.cutoff)) id.cutoff <- 95
      if(is.null(id.binsize)) id.binsize <- 0.1
   }else if(id.metric==2){
      if(ncol(rec)<6) stop("Requesting corrected identity, but .rec file doesn't have 6th column")
      id.reccol <- 6
      id.shortname <- 'cId.'
      id.fullname  <- 'Corrected identity'
      id.units     <- '%'
      id.hallmarks <- seq(0, 100, by=5)
      if(is.null(id.max)) id.max <- 100
      if(is.null(id.cutoff)) id.cutoff <- 95
      if(is.null(id.binsize)) id.binsize <- 0.1
   }else if(id.metric==3){
      id.reccol <- 4
      id.shortname <- 'BSc.'
      id.fullname  <- 'Bit score'
      id.units     <- 'bits'
      max.bs <- max(rec[, id.reccol])
      id.hallmarks <- seq(0, max.bs*1.2, by=50)
      if(is.null(id.max)) id.max <- max.bs
      if(is.null(id.cutoff)) id.cutoff <- 0.95 * max.bs
      if(is.null(id.binsize)) id.binsize <- 5
   }
   if(is.null(id.min)) id.min <- min(rec[, id.reccol]);
   if(is.null(pos.max)) pos.max <- max(lim[, 3]);
   id.lim <- c(id.min, id.max);
   pos.lim <- c(pos.min, pos.max)/1e6;
   id.breaks <- round((id.max-id.min)/id.binsize);
   pos.breaks <- round((pos.max-pos.min)/pos.binsize);
   if(is.null(main)) main <- paste('Recruitment plot of ', prefix, sep='');
   pos.marks=seq(pos.min, pos.max, length.out=pos.breaks+1)/1e6;
   id.marks=seq(id.min, id.max, length.out=id.breaks+1);
   id.topclasses <- 0;
   for(i in length(id.marks):1) if(id.marks[i]>id.cutoff) id.topclasses <- id.topclasses + 1;
   
   # Set-up image
   layout(matrix(c(3,4,1,2), nrow=2, byrow=TRUE), widths=c(2,1), heights=c(1,2));
   out <- list();

   # Recruitment plot
   if(verbose) cat("Rec. plot.\n")
   par(mar=c(5,4,0,0)+0.1);
   rec.hist <- matrix(0, nrow=pos.breaks, ncol=id.breaks);
   for(i in 1:nrow(rec)){
      id.class <- ceiling((id.breaks)*((rec[i, id.reccol]-id.min)/(id.max-id.min)));
      if(id.class<=id.breaks & id.class>0){
	 for(pos in rec[i, 1]:rec[i, 2]){
	    pos.class <- ceiling((pos.breaks)*((pos-pos.min)/(pos.max-pos.min)));
	    if(pos.class<=pos.breaks & pos.class>0) rec.hist[pos.class, id.class] <- rec.hist[pos.class, id.class]+1;
	 }
      }
   }
   id.top <- c((1-id.topclasses):0) + id.breaks;
   rec.col=gplots::colorpanel(256, rec.col1, rec.col2);
   image(x=pos.marks, y=id.marks, z=log10(rec.hist),
   		breaks=seq(0, log10(max(rec.hist)), length.out=1+length(rec.col)), col=rec.col,
		xlim=pos.lim, ylim=id.lim, xlab='Position in genome (Mbp)',
		ylab=paste(id.fullname, ' (',id.units,')', sep=''), xaxs='i', yaxs='r');
   if(!is.na(contig.col)) abline(v=c(lim$V2, lim$V3)/1e6, lty=1, col=contig.col);
   abline(h=id.hallmarks, lty=2, col=grey(0.7));
   abline(h=id.marks[id.top[1]], lty=3, col=grey(0.5))
   legend('bottomleft', 'Rec. plot', bg=rgb(1,1,1,2/3));
   out <- c(out, list(pos.marks=pos.marks, id.marks=id.marks));
   if(ret.recplot) out <- c(out, list(recplot=rec.hist));

   # Identity histogram
   if(verbose) cat(id.shortname, " hist.\n", sep='')
   par(mar=c(5,0,0,2)+0.1);
   id.hist <- id.summary.func(rec.hist);
   plot(1, t='n', xlim=c(1, max(id.hist)), ylim=id.lim, ylab='', yaxt='n', xlab=paste('Sequences (bp),', id.summary.name), log='x', ...);
   id.x <- rep(id.marks, each=2)[2:(id.breaks*2+1)]
   id.f <- rep(id.hist, each=2)[1:(id.breaks*2)]
   if(sum(id.f)>0){
      lines(id.f, id.x, lwd=ifelse(id.splines>0, 1/2, 2), type='o', pch='.');
      if(id.splines>0){
	 id.spline <- smooth.spline(id.x[id.f>0], log(id.f[id.f>0]), spar=id.splines)
	 lines(exp(id.spline$y), id.spline$x, lwd=2)
      }
   }
   
   abline(h=id.hallmarks, lty=2, col=grey(0.7));
   abline(h=id.marks[id.top[1]], lty=3, col=grey(0.5))
   legend('bottomright', paste(id.shortname, 'histogram'), bg=rgb(1,1,1,2/3));
   out <- c(out, list(id.mean=mean(rec[, id.reccol])));
   out <- c(out, list(id.median=median(rec[, id.reccol])));
   if(ret.mode)   out <- c(out, list(id.mode=modeest::mlv(rec[, id.reccol], method='mfv')$M));
   if(ret.hist)  out <- c(out, list(id.hist=id.hist));

   # Position histogram
   if(verbose) cat("Pos. hist.\n")
   par(mar=c(0,4,4,0)+0.1);
   h1<-rep(0,nrow(rec.hist)) ;
   h2<-rep(0,nrow(rec.hist)) ;
   pos.winsize <- (pos.max-pos.min+1)/pos.breaks;
   if(sum(rec.hist[, id.top])>0) h1 <- rowSums(matrix(rec.hist[, id.top], nrow=nrow(rec.hist)))/pos.winsize;
   if(sum(rec.hist[,-id.top])>0) h2 <- rowSums(matrix(rec.hist[,-id.top], nrow=nrow(rec.hist)))/pos.winsize;
   
   ymin <- min(1, h1[h1>0], h2[h2>0]);
   ymax <- max(10, h1, h2);
   if(is.na(ymin) || ymin<=0) ymin <- 1e-10;
   if(is.na(ymax) || ymax<=0) ymax <- 1;
   plot(1, t='n', xlab='', xaxt='n', ylab='Sequencing depth (X)', log='y', xlim=pos.lim,
   	ylim=c(ymin, ymax), xaxs='i', main=main, ...);
   if(!is.na(contig.col)) abline(v=c(lim[,2], lim[,3])/1e6, lty=1, col=contig.col);
   abline(h=10^c(0:5), lty=2, col=grey(0.7));
   if(sum(h2)>0){
      h2.x <- rep(pos.marks, each=2)[2:(pos.breaks*2+1)]
      h2.y <- rep(h2, each=2)[1:(pos.breaks*2)]
      lines(h2.x, h2.y, lwd=ifelse(pos.splines>0, 1/2, 2), col=grey(0.5));
      if(pos.splines>0){
         h2.spline <- smooth.spline(h2.x[h2.y>0], log(h2.y[h2.y>0]), spar=pos.splines)
	 lines(h2.spline$x, exp(h2.spline$y), lwd=2, col=grey(0.5))
      }
      if(ret.hist) out <- c(out, list(pos.hist.low=h2.y));
   }
   if(sum(h1)>0){
      h1.x <- rep(pos.marks, each=2)[2:(pos.breaks*2+1)]
      h1.y <- rep(h1, each=2)[1:(pos.breaks*2)]
      lines(h1.x, h1.y, lwd=ifelse(pos.splines>0, 1/2, 2), col=grey(0));
      if(pos.splines>0){
         h1.spline <- smooth.spline(h1.x[h1.y>0], log(h1.y[h1.y>0]), spar=pos.splines)
	 lines(h1.spline$x, exp(h1.spline$y), lwd=2, col=grey(0))
      }
      if(ret.hist) out <- c(out, list(pos.hist.top=h1.y));
   }
   legend('topleft', 'Pos. histogram', bg=rgb(1,1,1,2/3));
   out <- c(out, list(id.max=id.max, id.cutoff=id.marks[id.top[1]]));
   out <- c(out, list(seqdepth.mean.top=mean(h1)));
   out <- c(out, list(seqdepth.mean.low=mean(h2)));
   out <- c(out, list(seqdepth.mean=mean(h1+h2)));
   out <- c(out, list(seqdepth.median.top=median(h1)));
   out <- c(out, list(seqdepth.median.low=median(h2)));
   out <- c(out, list(seqdepth.median=median(h1+h2)));
   out <- c(out, list(id.metric=id.fullname));
   out <- c(out, list(id.summary=id.summary.name));
   
   # Legend
   par(mar=c(0,0,4,2)+0.1);
   plot(1, t='n', xlab='', xaxt='n', ylab='', yaxt='n', xlim=c(0,1), ylim=c(0,1), xaxs='r', yaxs='i', ...);
   text(1/2, 5/6, labels=paste('Reads per ', signif((pos.max-pos.min)/pos.breaks, 2), ' bp (rec. plot)', sep=''), pos=3);
   leg.col <- gplots::colorpanel(100, rec.col1, rec.col2);
   leg.lab <- signif(10^seq(0, log10(max(rec.hist)), length.out=10), 2);
   for(i in 1:10){
      for(j in 1:10){
         k <- (i-1)*10 + j;
	 polygon(c(k-1, k, k, k-1)/100, c(2/3, 2/3, 5/6, 5/6), border=leg.col[k], col=leg.col[k]);
      }
      text((i-0.5)/10, 2/3, labels=paste(leg.lab[i], ''), srt=90, pos=2, offset=0, cex=3/4);
   }
   legend('bottom',
   	legend=c('Contig boundary', 'Hallmark', paste(id.fullname, 'cutoff'),
		paste('Pos. hist.: ',id.shortname,' > ',signif(id.marks[id.top[1]],2),id.units,sep=''),
		paste('Pos. hist.: ',id.shortname,' < ',signif(id.marks[id.top[1]],2),id.units,sep='')), ncol=2,
   	col=grey(c(0.85, 0.7, 0.5, 0, 0.5)), lty=c(1,2,3,1,1), lwd=c(1,1,1,2,2), bty='n', inset=0.05, cex=5/6);
   return(out);
   ### A list with the following elements:
   ### 
   ### pos.marks: Midpoints of the position histogram.
   ### 
   ### id.matrix: Midpoints of the identity histogram.
   ### 
   ### recplot (if ret.recplot=TRUE): Matrix containing the recruitment plot values.
   ### 
   ### id.mean: Mean identity.
   ### 
   ### id.median: Median identity.
   ### 
   ### id.mode (if ret.mode=TRUE): Mode of the identity.
   ### 
   ### id.hist (if ret.hist=TRUE): Values of the identity histogram.
   ### 
   ### pos.hist.low (if ret.hist=TRUE): Values of the position histogram (depth) with "low"
   ### identity (i.e., below id.cutoff).
   ### 
   ### pos.hist.top (if ret.hist=TRUE): Values of the position histogram (depth) with "top"
   ### identity (i.e., above id.cutoff).
   ### 
   ### id.max: Value of id.max. This is returned because id.max=NULL may vary.
   ### 
   ### id.cutoff: Value of id.cutoff. This is returned because id.cutoff=NULL may vary.
   ### 
   ### seqdepth.mean.top: Average sequencing depth with identity above id.cutoff.
   ### 
   ### seqdepth.mean.low: Average sequencing depth with identity below id.cutoff.
   ### 
   ### seqdepth.mean.all: Average sequencing depth without identity filtering.
   ### 
   ### seqdepth.median.top: Median sequencing depth with identity above id.cutoff.
   ### 
   ### seqdepth.median.low: Median sequencing depth with identity below id.cutoff.
   ### 
   ### seqdepth.median.all: Median sequencing depth without identity filtering.
   ### 
   ### id.metric: Full name of the used identity metric.
   ### 
   ### id.summary: Full name of the summary method used to build the identity plot.
});

