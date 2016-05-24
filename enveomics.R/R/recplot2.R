#==============> Define S4 classes
setClass("enve.RecPlot2",
   ### Enve-omics representation of Recruitment plots. This object can
   ### be produced by `enve.recplot2` and supports S4 method plot.
   representation(
   counts='matrix',		##<< Counts as a two-dimensional histogram.
   pos.counts.in='numeric',	##<< Counts of in-group hits per position bin.
   pos.counts.out='numeric',	##<< Counts of out-group hits per position bin.
   id.counts='numeric',		##<< Counts per ID bin.
   id.breaks='numeric',		##<< Breaks of identity bins.
   pos.breaks='numeric',	##<< Breaks of position bins.
   seq.breaks='numeric',
   ### Limits of the subject sequences after concatenation.
   seq.names='character',	##<< Names of the subject sequences.
   id.metric='character',	##<< Metric used as 'identity'.
   id.ingroup='logical',	##<< Identity bins considered in-group.
   call='call')			##<< Call producing this object.
   ,package='enveomics.R'
   );
setClass("enve.RecPlot2.Peak",
### Enve-omics representation of a peak in the sequencing depth histogram
### of a Recruitment plot (see `enve.recplot2.findPeaks`).
   representation(
   dist='character',
   ### Distribution of the peak. Currently supported: 'norm' (normal) and 'sn' (skew-normal).
   values='numeric',
   ### Sequencing depth values predicted to conform the peak.
   values.res='numeric',
   ### Sequencing depth values not explained by this or previously identified peaks.
   mode='numeric',
   ### Seed-value of mode anchoring the peak.
   param.hat='list',
   ### Parameters of the distribution. A list of two values if dist='norm' (sd and mean),
   ### or three values if dist='sn' (omega=scale, alpha=shape, and xi=location). Note that the
   ### "dispersion" parameter is always first and the "location" parameter is always last.
   n.hat='numeric',
   ### Number of bins estimated to be explained by this peak. This should ideally be equal to
   ### the length of `values`, but it's not and integer.
   n.total='numeric',
   ### Total number of bins from which the peak was extracted.
   err.res='numeric',
   ### Error left after adding the peak.
   merge.logdist='numeric'
   ### Attempted `merge.logdist` parameter.
   ));
setMethod("$", "enve.RecPlot2", function(x, name) attr(x, name))
setMethod("$", "enve.RecPlot2.Peak", function(x, name) attr(x, name))

#==============> Define S4 methods
plot.enve.RecPlot2 <- function
   ### Plots an `enve.RecPlot2` object.
      (x,
      ### `enve.RecPlot2` object to plot.
      layout=matrix(c(5,5,2,1,4,3), nrow=2),
      ### Matrix indicating the position of the different panels in the layout, where:
      ### 0: Empty space, 1: Counts matrix, 2: position histogram (sequencing depth),
      ### 3: identity histogram, 4: Populations histogram (histogram of sequencing depths),
      ### 5: Color scale for the counts matrix (vertical), 6: Color scale of the counts
      ### matrix (horizontal). Only panels indicated here will be plotted. To plot only
      ### one panel simply set this to the number of the panel you want to plot.
      widths=c(1,7,2),
      ### Relative widths of the columns of `layout`.
      heights=c(1,2),
      ### Relative heights of the rows of `layout`.
      palette=grey((100:0)/100),
      ### Colors to be used to represent the counts matrix, sorted from no hits to the
      ### maximum sequencing depth.
      underlay.group=TRUE,
      ### If TRUE, it indicates the in-group and out-group areas couloured based on
      ### `in.col` and `out.col`. Requires support for semi-transparency.
      peaks.col='darkred',
      ### If not NA, it attempts to represent peaks in the population histogram in the
      ### specified color. Set to NA to avoid peak-finding.
      id.lim=range(x$id.breaks),
      ### Limits of identities to represent.
      pos.lim=range(x$pos.breaks),
      ### Limits of positions to represent (in bp, regardless of `pos.units`).
      pos.units=c('Mbp','Kbp','bp'),
      ### Units in which the positions should be represented (powers of 1,000 base pairs).
      mar=list('1'=c(5,4,1,1)+.1, '2'=c(ifelse(any(layout==1),1,5),4,4,1)+.1,
	    '3'=c(5,ifelse(any(layout==1),1,4),1,2)+0.1,
	    '4'=c(ifelse(any(layout==1),1,5),ifelse(any(layout==2),1,4),4,2)+0.1,
	    '5'=c(5,3,4,1)+0.1, '6'=c(5,4,4,2)+0.1),
      ### Margins of the panels as a list, with the character representation of the number
      ### of the panel as index (see `layout`).
      pos.splines=0,
      ### Smoothing parameter for the splines in the position histogram. Zero (0) for no splines.
      ### If non-zero, requires the stats package.
      id.splines=1/2,
      ### Smoothing parameter for the splines in the identity histogram. Zero (0) for no splines.
      ### If non-zero, requires the stats package.
      in.lwd=ifelse(pos.splines>0, 1/2, 2),
      ### Line width for the sequencing depth of in-group matches.
      out.lwd=ifelse(pos.splines>0, 1/2, 2),
      ### Line width for the sequencing depth of out-group matches.
      id.lwd=ifelse(id.splines>0, 1/2, 2),
      ### Line width for the identity histogram.
      in.col='darkblue',
      ### Color associated to in-group matches.
      out.col='lightblue',
      ### Color associated to out-group matches.
      id.col='black',
      ### Color for the identity histogram.
      breaks.col='#AAAAAA40',
      ### Color of the vertical lines indicating sequence breaks.
      peaks.opts=list(),
      ### Options passed to `enve.recplot2.findPeaks`, if `peaks.col` is not NA.
      ...
      ### Any other graphic parameters (currently ignored).
   ){
   pos.units	<- match.arg(pos.units);
   pos.factor	<- ifelse(pos.units=='bp',1,ifelse(pos.units=='Kbp',1e3,1e6));
   pos.lim	<- pos.lim/pos.factor;
   lmat <- layout;
   for(i in 1:6) if(!any(layout==i)) lmat[layout>i] <- lmat[layout>i]-1;

   layout(lmat, widths=widths, heights=heights);
   ori.mar <- par('mar');

   # Essential vars
   counts	<- x$counts;

   id.ingroup	<- x$id.ingroup;
   id.counts	<- x$id.counts;
   id.breaks	<- x$id.breaks;
   id.mids	<- (id.breaks[-length(id.breaks)]+id.breaks[-1])/2;
   id.binsize	<- id.breaks[-1] - id.breaks[-length(id.breaks)];

   pos.counts.in	<- x$pos.counts.in;
   pos.counts.out	<- x$pos.counts.out;
   pos.breaks		<- x$pos.breaks/pos.factor;
   pos.mids		<- (pos.breaks[-length(pos.breaks)]+pos.breaks[-1])/2;
   pos.binsize		<- (pos.breaks[-1] - pos.breaks[-length(pos.breaks)])*pos.factor;

   seqdepth.in  <- pos.counts.in/pos.binsize;
   seqdepth.out <- pos.counts.out/pos.binsize;
   seqdepth.lim <- range(c(seqdepth.in[seqdepth.in>0],seqdepth.out[seqdepth.out>0]))*c(1/2,2);

   if(underlay.group){
      in.bg  <- do.call(rgb, c(as.list(col2rgb(in.col)), list(maxColorValue=256, alpha=62)));
      out.bg <- do.call(rgb, c(as.list(col2rgb(out.col)[,1]), list(maxColorValue=256, alpha=52)));
   }

   # Counts matrix
   if(any(layout==1)){
      par(mar=mar[['1']]);
      plot(1, t='n', bty='l',
	 xlim=pos.lim, xlab=paste('Position in genome (',pos.units,')',sep=''), xaxs='i',
	    ylim=id.lim,  ylab=x$id.metric, yaxs='i');
      if(underlay.group){
	 rect(pos.lim[1], id.lim[1], pos.lim[2], min(id.breaks[c(id.ingroup,TRUE)]), col=out.bg, border=NA);
	 rect(pos.lim[1], min(id.breaks[c(id.ingroup,TRUE)]), pos.lim[2], id.lim[2], col=in.bg,  border=NA);
      }
      abline(v=x$seq.breaks/pos.factor, col=breaks.col);
      image(x=pos.breaks, y=id.breaks, z=log10(counts),col=palette, bg=grey(1,0),
	 breaks=seq(-.1,log10(max(counts)), length.out=1+length(palette)), add=TRUE);
   }

   # Position histogram
   if(any(layout==2)){
      par(mar=mar[['2']]);
      if(any(layout==1)){
	 xlab=''
	 xaxt='n'
      }else{
	 xlab=paste('Position in genome (',pos.units,')',sep='')
	 xaxt='s'
      }
      plot(1,t='n', bty='l', log='y',
	 xlim=pos.lim, xlab=xlab, xaxt=xaxt, xaxs='i',
	 ylim=seqdepth.lim, yaxs='i', ylab='Sequencing depth (X)');
      abline(v=x$seq.breaks/pos.factor, col=breaks.col)
      pos.x <- rep(pos.breaks,each=2)[-c(1,2*length(pos.breaks))]
      pos.f <- rep(seqdepth.in,each=2)
      lines(pos.x, rep(seqdepth.out,each=2), lwd=out.lwd, col=out.col);
      lines(pos.x, pos.f, lwd=in.lwd, col=in.col);
      if(pos.splines > 0){
	 pos.spline <- smooth.spline(pos.x[pos.f>0], log(pos.f[pos.f>0]), spar=pos.splines)
	 lines(pos.spline$x, exp(pos.spline$y), lwd=2, col=in.col)
      }
      if(any(pos.counts.out==0)) rect(pos.breaks[c(pos.counts.out==0,FALSE)], seqdepth.lim[1], pos.breaks[c(FALSE,pos.counts.out==0)], seqdepth.lim[1]*3/2, col=out.col, border=NA);
      if(any(pos.counts.in==0))  rect(pos.breaks[c(pos.counts.in==0,FALSE)],  seqdepth.lim[1], pos.breaks[c(FALSE,pos.counts.in==0)],  seqdepth.lim[1]*3/2, col=in.col,  border=NA);
   }

   # Identity histogram
   if(any(layout==3)){
      par(mar=mar[['3']]);
      if(any(layout==1)){
	 ylab=''
	 yaxt='n'
      }else{
	 ylab=x$id.metric
	 yaxt='s'
      }
      if(sum(id.counts>0) >= 4){
	 id.counts.range <- range(id.counts[id.counts>0])*c(1/2,2);
	 plot(1,t='n', bty='l', log='x',
	       xlim=id.counts.range, xlab='bps per bin', xaxs='i',
	       ylim=id.lim, yaxs='i', ylab=ylab, yaxt=yaxt);
	 if(underlay.group){
	    rect(id.counts.range[1], id.lim[1], id.counts.range[2], min(id.breaks[c(id.ingroup,TRUE)]), col=out.bg, border=NA);
	    rect(id.counts.range[1], min(id.breaks[c(id.ingroup,TRUE)]), id.counts.range[2], id.lim[2], col=in.bg,  border=NA);
	 }
	 id.f <- rep(id.counts,each=2)
	 id.x <- rep(id.breaks,each=2)[-c(1,2*length(id.breaks))]
	 lines(id.f, id.x, lwd=id.lwd, col=id.col);
	 if(id.splines > 0){
	    id.spline <- smooth.spline(id.x[id.f>0], log(id.f[id.f>0]), spar=id.splines)
	    lines(exp(id.spline$y), id.spline$x, lwd=2, col=id.col)
	 }
      }else{
	 plot(1,t='n',bty='l',xlab='', xaxt='n', ylab='', yaxt='n')
	 text(1,1,labels='Insufficient data', srt=90)
      }
   }

   # Populations histogram
   peaks <- NA;
   if(any(layout==4)){
      par(mar=mar[['4']]);
      if(any(layout==2)){
	 ylab=''
	 yaxt='n'
      }else{
	 ylab='Sequencing depth (X)'
	 yaxt='s'
      }
      h.breaks <- seq(log10(seqdepth.lim[1]*2), log10(seqdepth.lim[2]/2), length.out=200);
      h.in <- hist(log10(seqdepth.in), breaks=h.breaks, plot=FALSE);
      h.out <- hist(log10(seqdepth.out), breaks=h.breaks, plot=FALSE);
      plot(1, t='n', log='y',
	 xlim=range(c(h.in$counts,h.out$counts,sum(pos.counts.in==0))), xaxs='r', xlab='', xaxt='n',
	 ylim=seqdepth.lim, yaxs='i', ylab=ylab, yaxt=yaxt);
      y.tmp.in <- c(rep(10^h.in$breaks,each=2),seqdepth.lim[1]*c(1,1,3/2,3/2));
      y.tmp.out <- c(rep(10^h.out$breaks,each=2),seqdepth.lim[1]*c(1,1,3/2,3/2));
      lines(c(0,rep(h.out$counts,each=2),0,0,rep(sum(pos.counts.out==0),2),0), y.tmp.out, col=out.col);
      polygon(c(0,rep(h.in$counts,each=2),0,0,rep(sum(pos.counts.in==0),2),0), y.tmp.in, border=NA, col=in.col);
      if(!is.na(peaks.col)){
	 o	<- peaks.opts; o$x = x;
	 peaks	<- do.call(enve.recplot2.findPeaks, o);
	 h.mids <- (10^h.breaks[-1] + 10^h.breaks[-length(h.breaks)])/2
	 if(!is.null(peaks) & length(peaks)>0){
	    pf <- h.mids*0;
	    for(i in 1:length(peaks)){
	       cnt <- enve.recplot2.__peakHist(peaks[[i]], h.mids)
	       lines(cnt, h.mids, col='red');
	       pf <- pf+cnt;
	       axis(4, at=peaks[[i]]$param.hat[[ length(peaks[[i]]$param.hat) ]], letters[i], las=1, hadj=1/2)
	    }
	    lines(pf, h.mids, col='red',lwd=1.5);
	    legend('bottomright', legend=paste(
	       letters[1:length(peaks)],'. ',
	       signif(as.numeric(lapply(peaks, function(x) tail(as.numeric(x$param.hat),n=1))),3),'X (',
	       signif(100*as.numeric(lapply(peaks, function(x) (length(x$values)/x$n.total))), 3), '%, err: ',
	       signif(as.numeric(lapply(peaks, function(x) x$err.res)), 3), ')',
	       sep=''), bty='n');
	 }
      }
   }

   # Color scale
   count.bins <- 10^seq(log10(min(counts[counts>0])), log10(max(counts)), length.out=1+length(palette));
   if(any(layout==5)){
      par(mar=mar[['5']]);
      plot(1,t='n',log='y',xlim=0:1,xaxt='n',xlab='',xaxs='i',ylim=range(count.bins), yaxs='i', ylab='');
      rect(0,count.bins[-length(count.bins)],1,count.bins[-1],col=palette,border=NA);
   }
   if(any(layout==6)){
      par(mar=mar[['6']]);
      plot(1,t='n',log='x',ylim=0:1,yaxt='n',ylab='',yaxs='i',xlim=range(count.bins), xaxs='i',xlab='');
      rect(count.bins[-length(count.bins)],0,count.bins[-1],1,col=palette,border=NA);
   }
   
   par(mar=ori.mar);
   return(peaks);
   ### Returns a list of `enve.RecPlot2.Peak` objects (see `enve.recplot2.findPeaks`). If `peaks.col`=NA or
   ### `layout` doesn't include 4, returns NA.
}

#==============> Define core functions
enve.recplot2 <- function(
   ### Produces recruitment plots provided that BlastTab.catsbj.pl has
   ### been previously executed.
      prefix,
      ### Path to the prefix of the BlastTab.catsbj.pl output files. At
      ### least the files .rec and .lim must exist with this prefix.
      plot=TRUE,
      ### Should the object be plotted?
      pos.breaks=1e3,
      ### Breaks in the positions histogram. It can also be a vector of break
      ### points, and values outside the range are ignored. If zero (0), it
      ### uses the sequence breaks as defined in the .lim file, which means
      ### one bin per contig (or gene, if the mapping is agains genes).
      id.breaks=300,
      ### Breaks in the identity histogram. It can also be a vector of break
      ### points, and values outside the range are ignored.
      id.metric=c('identity', 'corrected identity', 'bit score'),
      ### Metric of identity to be used (Y-axis). Corrected identity is only
      ### supported if the original BLAST file included sequence lengths.
      id.summary=sum,
      ### Function summarizing the identity bins. Other recommended options
      ### include: `median` to estimate the median instead of total bins, and
      ### `function(x) mlv(x,method='parzen')$M` to estimate the mode.
      id.cutoff=95,
      ### Cutoff of identity metric above which the hits are considered
      ### 'in-group'. The 95% identity corresponds to the expectation of
      ### ANI<95% within species.
      threads=2,
      ### Number of threads to use.
      verbose=TRUE,
      ### Indicates if the function should report the advance.
      ...
      ### Any additional parameters supported by `plot.enve.RecPlot2`.
   ){
   # Settings
   id.metric <- match.arg(id.metric);
   
   #Read files
   if(verbose) cat("Reading files.\n")
   rec <- read.table(paste(prefix, ".rec", sep=""), sep="\t", comment.char="",
      quote="");
   lim <- read.table(paste(prefix, ".lim", sep=""), sep="\t", comment.char="",
      quote="", as.is=TRUE);
   
   # Build matrix
   if(verbose) cat("Building counts matrix.\n")
   if(id.metric=="corrected identity" & ncol(rec)<6){
      stop("Requesting corrected identity, but .rec file doesn't have 6th column")
   }
   rec.idcol <- ifelse(id.metric=="identity", 3, ifelse(id.metric=="corrected identity", 6, 4));
   if(length(pos.breaks)==1){
      if(pos.breaks>0){
         pos.breaks <- seq(min(lim[,2]), max(lim[,3]), length.out=pos.breaks+1);
      }else{
         pos.breaks <- c(lim[,2], tail(lim[,3], n=1))
      }
   }
   if(length(id.breaks)==1){
      id.breaks <- seq(min(rec[,rec.idcol]), max(rec[,rec.idcol]),
	 length.out=id.breaks+1);
   }
   
   # Run in parallel
   cl		<- makeCluster(threads)
   rec.l	<- list()
   thl		<- ceiling(nrow(rec)/threads)
   for(i in 0:(threads-1)){
      rec.l[[i+1]] <- list(rec=rec[ (i*thl+1):min(((i+1)*thl),nrow(rec)), ],
			verbose=ifelse(i==0, verbose, FALSE))
   }
   counts.l	<- clusterApply(cl, rec.l, enve.recplot2.__counts,
			pos.breaks=pos.breaks, id.breaks=id.breaks,
			rec.idcol=rec.idcol)
   counts	<- counts.l[[1]]
   if(threads>1) for(i in 2:threads) counts <- counts + counts.l[[i]]
   stopCluster(cl)
   
   # Estimate 1D histograms
   if(verbose) cat("Building histograms.\n")
   id.mids	<- (id.breaks[-length(id.breaks)]+id.breaks[-1])/2;
   id.ingroup	<- (id.mids > id.cutoff);
   id.counts	<- apply(counts, 2, id.summary);
   pos.counts.in   <- apply(counts[,id.ingroup], 1, sum);
   pos.counts.out  <- apply(counts[,!id.ingroup], 1, sum);

   # Plot and return
   recplot <- new('enve.RecPlot2',
      counts=counts, id.counts=id.counts, pos.counts.in=pos.counts.in,
      pos.counts.out=pos.counts.out,
      id.breaks=id.breaks, pos.breaks=pos.breaks,
      seq.breaks=c(lim[1,2], lim[,3]), seq.names=lim[,1],
      id.ingroup=id.ingroup,id.metric=id.metric,
      call=match.call());
   if(plot){
      if(verbose) cat("Plotting.\n")
      plot(recplot, ...);
   }
   return(recplot);
   ### Returns an object of class `enve.RecPlot2`.
}

enve.recplot2.findPeaks <- function(
   ### Identifies peaks in the population histogram potentially indicating
   ### sub-population mixtures.
      x,
      ### An `enve.RecPlot2` object.
      min.points=50,
      ### Minimum number of points in the quantile-estimation-range (`quant.est`) to estimate a peak.
      quant.est=c(0.005, 0.995),
      ### Range of quantiles to be used in the estimation of a peak's parameters.
      mlv.opts=list(method='parzen'),
      ### Options passed to `mlv` to estimate the mode.
      fitdist.opts.sn=list(distr='sn', method='qme', probs=c(.1, .5, .8), start=list(omega=1, alpha=-1), lower=c(1e-6, -Inf, 0), upper=c(Inf, 0, Inf)),
      ### Options passed to `fitdist` to estimate the standard deviation if with.skewness=TRUE. Note that
      ### the `start` parameter will be ammended with xi=estimated mode for each peak.
      fitdist.opts.norm=list(distr='norm', method='qme', probs=c(.4,.6), start=list(sd=1), lower=c(1e-8, 0)),
      ### Options passed to `fitdist` to estimate the standard deviation if with.skewness=FALSE. Note that
      ### the `start` parameter will be ammended with mean=estimated mode for each peak.
      rm.top=0.05,
      ### Top-values to remove before finding peaks, as a quantile probability. This step
      ### is useful to remove highly conserved regions, but can be turned off by setting
      ### rm.top=0. The quantile is determined *after* removing zero-coverage windows.
      with.skewness=TRUE,
      ### Allow skewness correction of the peaks. Typically, the sequencing-depth distribution
      ### for a single peak is left-skewed, due partly (but not exclusively) to fragmentation
      ### and mapping sensitivity. See Lindner et al 2013, Bioinformatics 29(10):1260-7 for
      ### an alternative solution for the first problem (fragmentation) called "tail
      ### distribution".
      optim.rounds=200,
      ### Maximum rounds of peak optimization.
      optim.epsilon=1e-8,
      ### Trace change at which optimization stops (unless `optim.rounds` is reached first). The
      ### trace change is estimated as the sum of square differences between parameters in one round and
      ### those from two rounds earlier (to avoid infinite loops from approximation).
      merge.logdist=log(1.25),
      ### Maximum value of |log-ratio| between centrality parameters in peaks to attempt merging. The default
      ### of ~0.22 corresponds to a maximum difference of 25%.
      verbose=FALSE
      ### Display (mostly debugging) information.
   ){
   
   # Essential vars
   pos.binsize	<- x$pos.breaks[-1] - x$pos.breaks[-length(x$pos.breaks)];
   seqdepth.in	<- x$pos.counts.in/pos.binsize;
   lsd1 <- seqdepth.in[seqdepth.in>0];
   lsd1 <- lsd1[ lsd1 < quantile(lsd1, 1-rm.top, names=FALSE) ]
   if(with.skewness){ fitdist.opts <- fitdist.opts.sn }else{ fitdist.opts <- fitdist.opts.norm }
   peaks.opts <- list(lsd1=lsd1, min.points=min.points, quant.est=quant.est, mlv.opts=mlv.opts,
      fitdist.opts=fitdist.opts, with.skewness=with.skewness, optim.rounds=optim.rounds,
      optim.epsilon=optim.epsilon, verbose=verbose, n.total=length(lsd1), merge.logdist=merge.logdist)
   
   # Find seed peaks
   if(verbose) cat('Mowing peaks for n =',length(lsd1),'\n')
   peaks <- enve.recplot2.__findPeaks(peaks.opts);

   # Merge overlapping peaks
   if(verbose) cat('Trying to merge',length(peaks),'peaks\n')
   merged <- (length(peaks)>1)
   while(merged){
      merged <- FALSE
      ignore <- c()
      peaks2 <- list();
      for(i in 1:length(peaks)){
	 if(i %in% ignore) next
	 p <- peaks[[ i ]]
	 j <- enve.recplot2.__whichClosestPeak(p, peaks)
	 p2 <- peaks[[ j ]]
	 if( abs(log(p$param.hat[[ length(p$param.hat) ]]/p2$param.hat[[ length(p2$param.hat) ]])) < merge.logdist ){
	    if(verbose) cat('==> Attempting a merge at',p$param.hat[[ length(p$param.hat) ]],'&',p2$param.hat[[ length(p2$param.hat) ]],'X\n');
	    peaks.opts$lsd1 <- c(p$values, p2$values)
	    p.new <- enve.recplot2.__findPeaks(peaks.opts)
	    if(length(p.new)==1){
	       peaks2[[ length(peaks2)+1 ]] <- p.new[[ 1 ]]
	       ignore <- c(ignore, j)
	       merged <- TRUE
	    }
	 }
	 if(!merged) peaks2[[ length(peaks2)+1 ]] <- p
      }
      peaks <- peaks2
      if(length(peaks)==1) break
   }
   
   if(verbose) cat('Found',length(peaks),'peak(s)\n')
   return(peaks);
   ### Returns a list of `enve.RecPlot2.Peak` objects.
}

enve.recplot2.corePeak <- function
   ### Finds the peak in a list of peaks that is most likely to represent the
   ### "core genome" of a population.
      (x
      ### `list` of `enve.RecPlot2.Peak` objects.
   ){
   # Find the peak with maximum depth (centrality)
   maxPeak <- x[[
	 which.max(as.numeric(lapply(x,
	    function(y) y$param.hat[[ length(y$param.hat) ]])))
      ]]
   # If a "larger" peak (a peak explaining more bins of the genome) is within
   # the "merge.logdist" distance, take that one instead.
   corePeak <- maxPeak
   for(p in x){
      sz.d = log(length(p$values)/length(corePeak$values))
      if(sz.d < 0)
	 next;
      sq.d = log(p$param.hat[[ length(p$param.hat) ]]/maxPeak$param.hat[[ length(maxPeak$param.hat) ]] )
      if(abs(sq.d) < maxPeak$merge.logdist+sz.d/5)
         corePeak <- p
   }
   return(corePeak)
}


enve.recplot2.__counts <- function
   ### Internal ancilliary function (see `enve.recplot2`).
      (x, pos.breaks, id.breaks, rec.idcol){
   rec <- x$rec
   verbose <- x$verbose
   counts <- matrix(0, nrow=length(pos.breaks)-1, ncol=length(id.breaks)-1);
   for(i in 1:nrow(rec)){
      if(verbose & i%%100==0) cat("   [",signif(i*100/nrow(rec),3),"% ]   \r");
      y.bin <- which(rec[i,rec.idcol]>=id.breaks[-length(id.breaks)] & rec[i,rec.idcol]<=id.breaks[-1])[1] ;
      for(pos in rec[i,1]:rec[i,2]){
	 x.bin <- which(pos>=pos.breaks[-length(pos.breaks)] & pos<=pos.breaks[-1])[1] ;
	 counts[x.bin, y.bin] <- counts[x.bin, y.bin]+1 ;
      }
   }
   return(counts);
}

enve.recplot2.__peakHist <- function
   ### Internal ancilliary function (see `enve.RecPlot2.Peak`).
      (x, mids, counts=TRUE){
   d.o <- x$param.hat
   d.o$x <- mids
   prob  <- do.call(paste('d', x$dist, sep=''), d.o)
   if(!counts) return(prob)
   if(length(x$values)>0) return(prob*length(x$values)/sum(prob))
   return(prob*x$n.hat/sum(prob))
}

enve.recplot2.__findPeak <- function
   ### Internall ancilliary function (see `enve.recplot2.findPeaks`).
      (lsd1, min.points, quant.est, mlv.opts, fitdist.opts, with.skewness,
      optim.rounds, optim.epsilon, n.total, merge.logdist, verbose
   ){
   dist	<- ifelse(with.skewness, 'sn', 'norm');
   
   # Find peak
   o <- mlv.opts; o$x = lsd1;
   mode1 <- do.call(mlv, o)$M;
   if(verbose) cat('Anchoring at mode =',mode1,'\n')
   param.hat <- fitdist.opts$start; last.hat <- param.hat;
   lim <- NA;
   if(with.skewness){ param.hat$xi <- mode1 }else{ param.hat$mean <- mode1 }
   
   # Refine peak parameters
   for(round in 1:optim.rounds){
      param.hat[[ 1 ]] <- param.hat[[ 1 ]]/diff(quant.est) # <- expand dispersion
      lim.o <- param.hat
      lim.o$p <- quant.est; lim <- do.call(paste('q',dist,sep=''), lim.o)
      lsd1.pop <- lsd1[(lsd1>lim[1]) & (lsd1<lim[2])];
      if(verbose) cat(' Round', round, 'with n =',length(lsd1.pop),'and params =',as.numeric(param.hat),' \r')
      if(length(lsd1.pop) < min.points) break;
      o <- fitdist.opts; o$data = lsd1.pop; o$start = param.hat;
      last.last.hat <- last.hat
      last.hat <- param.hat
      param.hat <- as.list(do.call(fitdist, o)$estimate);
      if(any(is.na(param.hat))){
	 if(round>1) param.hat <- last.hat;
	 break;
      }
      if(round>2) if(sum((as.numeric(last.last.hat)-as.numeric(param.hat))^2) < optim.epsilon) break;
   }
   if(verbose) cat('\n')
   if(is.na(param.hat[1]) | is.na(lim[1])) return(NULL);

   # Mow distribution
   lsd2 <- c();
   lsd.pop <- c();
   n.hat <- length(lsd1.pop)/diff(quant.est)
   peak <- new('enve.RecPlot2.Peak', dist=dist, values=as.numeric(), mode=mode1, param.hat=param.hat,
      n.hat=n.hat, n.total=n.total, merge.logdist=merge.logdist)
   peak.breaks <- seq(min(lsd1), max(lsd1), length=20)
   peak.cnt <- enve.recplot2.__peakHist(peak, (peak.breaks[-length(peak.breaks)]+peak.breaks[-1])/2)
   for(i in 2:length(peak.breaks)){
      values <- lsd1[ (lsd1 >= peak.breaks[i-1]) & (lsd1 < peak.breaks[i]) ]
      n.exp <- peak.cnt[i-1]
      if(n.exp==0) n.exp=0.1
      if(length(values)==0) next
      in.peak <- runif(length(values)) <= n.exp/length(values)
      lsd2 <- c(lsd2, values[!in.peak])
      lsd.pop <- c(lsd.pop, values[in.peak])
   }
   if(length(lsd.pop) < min.points) return(NULL)

   # Return peak
   attr(peak, 'values') <- lsd.pop
   attr(peak, 'values.res') <- lsd2
   attr(peak, 'err.res') <- 1-(cor(hist(lsd.pop, breaks=peak.breaks, plot=FALSE)$counts, hist(lsd1, breaks=peak.breaks, plot=FALSE)$counts)+1)/2
   if(verbose) cat(' Extracted peak with n =',length(lsd.pop),'with expected n =',n.hat,'\n')
   return(peak)
}

enve.recplot2.__findPeaks <- function
   ### Internal ancilliary function (see `enve.recplot2.findPeaks`).
      (peaks.opts){
   peaks <- list()
   while(length(peaks.opts$lsd1) > peaks.opts$min.points){
      peak <- do.call(enve.recplot2.__findPeak, peaks.opts)
      if(is.null(peak)) break
      peaks[[ length(peaks)+1 ]] <- peak
      peaks.opts$lsd1 <- peak$values.res
   }
   return(peaks)
}


enve.recplot2.__whichClosestPeak <- function
   ### Internal ancilliary function (see `enve.recplot2.findPeaks`).
      (peak, peaks){
   dist <- as.numeric(lapply(peaks, function(x) abs(log(x$param.hat[[ length(x$param.hat) ]]/peak$param.hat[[ length(peak$param.hat) ]] ))))
   dist[ dist==0 ] <- Inf
   return(which.min(dist))
}

