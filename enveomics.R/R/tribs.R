
# Use as:
# > # Estimate reference (null) model:
# > tab <- read.table('Ecoli-ML-dmatrix.txt', sep='\t', h=T, row.names=1)
# > dist <- as.dist(tab);
# > all.dist <- enve.tribs(dist);
# > 
# > # Estimate subset (test) model:
# > lee <- read.table('LEE-strains.txt', as.is=T)$V1
# > lee.dist <- enve.tribs(dist, lee, subsamples=seq(0,1,by=0.05), threads=12,
# +    verbosity=2, pre.tribs=all.dist.merge);
# ...
# > 
# > # Plot reference and selection at different subsampling levels:
# > plot(all.dist, t='boxplot');
# > plot(lee, new=FALSE, col='darkred');
# ...
# > 
# > # Test significance of overclustering (or overdispersion):
# > lee.test <- enve.tribs.test(dist, lee, pre.tribs=all.dist.merge,
# +    verbosity=2, threads=12);
# > summary(lee.test);
# > plot(lee.test);
# ...


#==============> Define S4 classes
setClass("enve.TRIBS",
   ### Enve-omics representation of "Transformed-space Resampling In Biased Sets
   ### (TRIBS)".  This object represents sets of distances between objects,
   ### sampled nearly-uniformly at random in "distance space".  Subsampling
   ### without selection is trivial, since both the distances space and the
   ### selection occur in the same transformed space. However, it's useful to
   ### compare randomly subsampled sets against a selected set of objects. This
   ### is intended to identify overdispersion or overclustering (see
   ### `enve.TRIBStest`) of a subset against the entire collection of objects
   ### with minimum impact of sampling biases. This object can be produced by
   ### `enve.tribs` and supports S4 methods `plot` and `summary`.
   representation(
   distance='numeric',	##<< Centrality measurement of the distances between the
			##<< selected objects (without subsampling).
   points='matrix',	##<< Position of the different objects in distance
			##<< space.
   distances='matrix',	##<< Subsampled distances, where the rows are replicates
			##<< and the columns are subsampling levels.
   spaceSize='numeric',	##<< Number of objects.
   selSize='numeric',	##<< Number of selected objects.
   dimensions='numeric',##<< Number of dimensions in the distance space.
   subsamples='numeric',##<< Subsampling levels (as fractions, from 0 to 1).
   call='call')		##<< Call producing this object.
   ,package='enveomics.R'
   );
setClass("enve.TRIBStest",
   ### Test of significance of overclustering or overdispersion in a selected
   ### set of objects with respect to the entire set (see `enve.TRIBS`). This
   ### object can be produced by `enve.tribs.test` and supports S4 methods
   ### `plot` and `summary`.
   representation(
   pval.gt='numeric',
   ### P-value for the overdispersion test.
   pval.lt='numeric',
   ### P-value for the overclustering test.
   all.dist='numeric',
   ### Empiric PDF of distances for the entire dataset (subsampled at selection
   ### size).
   sel.dist='numeric',
   ### Empiric PDF of distances for the selected objects (without subsampling).
   diff.dist='numeric',
   ### Empiric PDF of the difference between `all.dist` and `sel.dist`. The
   ### p-values are estimating by comparing areas in this PDF greater than and
   ### lesser than zero.
   dist.mids='numeric',
   ### Midpoints of the empiric PDFs of distances.
   diff.mids='numeric',
   ### Midpoints of the empiric PDF of difference of distances.
   call='call')
   ### Call producing this object.
   ,package='enveomics.R'
   );

#==============> Define S4 methods
summary.enve.TRIBS <- function
   ### Summary of an `enve.TRIBS` object.
   (object,
   ### `enve.TRIBS` object.
   ...
   ### No additional parameters are currently supported.
   ){
   cat('===[ enve.TRIBS ]-------------------------\n');
   cat('Selected',attr(object,'selSize'),'of',
      attr(object,'spaceSize'),'objects in',
      attr(object,'dimensions'),'dimensions.\n');
   cat('Collected',length(attr(object,'subsamples')),'subsamples with',
      nrow(attr(object,'distances')),'replicates each.\n');
   cat('------------------------------------------\n');
   cat('call:',as.character(attr(object,'call')),'\n');
   cat('------------------------------------------\n');
}

plot.enve.TRIBS <- function
   ### Plot an `enve.TRIBS` object.
   (x,
   ### `enve.TRIBS` object to plot.
   new=TRUE,
   ### Should a new canvas be drawn?
   type=c('boxplot', 'points'),
   ### Type of plot. The 'points' plot shows all the replicates, the 'boxplot'
   ### plot represents the values found by `boxplot.stats` as areas, and plots
   ### the outliers as points.
   col='#00000044',
   ### Color of the areas and/or the points.
   pt.cex=1/2,
   ### Size of the points.
   pt.pch=19,
   ### Points character.
   pt.col=col,
   ### Color of the points.
   ln.col=col,
   ### Color of the lines.
   ...
   ### Any additional parameters supported by `plot`.
   ){
   type <- match.arg(type);
   plot.opts <- list(xlim=range(attr(x,'subsamples'))*attr(x,'selSize'),
      ylim=range(attr(x,'distances')), ..., t='n', x=1);
   if(new) do.call(plot, plot.opts);
   abline(h=attr(x,'distance'), lty=3, col=ln.col);
   replicates <- nrow(attr(x,'distances'));
   if(type=='points'){
      for(i in 1:ncol(attr(x,'distances')))
	 points(rep(round(attr(x,'subsamples')[i]*attr(x,'selSize')),
	    replicates), attr(x,'distances')[,i], cex=pt.cex, pch=pt.pch,
	    col=pt.col);
   }else{
      stats <- matrix(NA, nrow=7, ncol=ncol(attr(x,'distances')));
      for(i in 1:ncol(attr(x,'distances'))){
	 b <- boxplot.stats(attr(x,'distances')[,i]);
	 points(rep(round(attr(x,'subsamples')[i]*attr(x,'selSize')),
	    length(b$out)), b$out, cex=pt.cex, pch=pt.pch, col=pt.col);
	 stats[, i] <- c(b$conf, b$stats[c(1,5,2,4,3)]);
      }
      x <- round(attr(x,'subsamples')*attr(x,'selSize'))
      for(i in c(1,3,5))
	 polygon(c(x, rev(x)), c(stats[i,], rev(stats[i+1,])), border=NA,
	    col=col);
      lines(x, stats[7,], col=ln.col, lwd=2);
   }
}

summary.enve.TRIBStest <- function
   ### Summary of an `enve.TRIBStest` object.
   (object,
   ### `enve.TRIBStest` object.
   ...
   ### No additional parameters are currently supported.
   ){
   cat('===[ enve.TRIBStest ]---------------------\n');
   cat('Alternative hypothesis:\n');
   cat('   The distances in the selection are\n');
   if(attr(object, 'pval.gt') > attr(object, 'pval.lt')){
      cat('   smaller than in the entire dataset\n   (overclustering)\n');
   }else{
      cat('   larger than in the entire dataset\n   (overdispersion)\n');
   }
   p.val <- min(attr(object, 'pval.gt'), attr(object, 'pval.lt'));
   if(p.val==0){
      diff.dist <- attr(object, 'diff.dist');
      p.val.lim <- min(diff.dist[diff.dist>0]);
      cat('\n   P-value <= ', signif(p.val.lim, 4), sep='');
   }else{
      p.val.lim <- p.val;
      cat('\n   P-value: ', signif(p.val, 4), sep='');
   }
   cat(' ', ifelse(p.val.lim<=0.01, "**", ifelse(p.val.lim<=0.05, "*", "")),
      '\n', sep='');
   cat('------------------------------------------\n');
   cat('call:',as.character(attr(object,'call')),'\n');
   cat('------------------------------------------\n');
}

plot.enve.TRIBStest <- function
   ### Plots an `enve.TRIBStest` object.
   (x,
   ### `enve.TRIBStest` object to plot.
   type=c('overlap', 'difference'),
   ### What to plot. 'overlap' generates a plot of the two contrasting empirical
   ### PDFs (to compare against each other), 'difference' produces a plot of the
   ### differences between the empirical PDFs (to compare against zero).
   col='#00000044',
   ### Main color of the plot if type='difference'.
   col1=col,
   ### First color of the plot if type='overlap'.
   col2='#44001144',
   ### Second color of the plot if type='overlap'.
   ylab='Probability',
   ### Y-axis label.
   xlim=range(attr(x, 'dist.mids')),
   ### X-axis limits.
   ylim=c(0,max(c(attr(x, 'all.dist'), attr(x, 'sel.dist')))),
   ### Y-axis limits.
   ...
   ### Any other graphical arguments.
   ){
   type <- match.arg(type);
   if(type=='overlap'){
      plot.opts <- list(xlim=xlim, ylim=ylim, ylab=ylab, ..., t='n', x=1);
      do.call(plot, plot.opts);
      bins <- length(attr(x, 'dist.mids'))
      polygon(attr(x, 'dist.mids')[c(1, 1:bins, bins)],
	 c(0,attr(x, 'all.dist'),0), col=col1,
	 border=do.call(rgb, as.list(c(col2rgb(col1)/256, 0.5))));
      polygon(attr(x, 'dist.mids')[c(1, 1:bins, bins)],
	 c(0,attr(x, 'sel.dist'),0), col=col2,
	 border=do.call(rgb, as.list(c(col2rgb(col2)/256, 0.5))));
   }else{
      plot.opts <- list(xlim=range(attr(x, 'diff.mids')),
	 ylim=c(0,max(attr(x, 'diff.dist'))), ylab=ylab, ..., t='n', x=1);
      do.call(plot, plot.opts);
      bins <- length(attr(x, 'diff.mids'));
      polygon(attr(x, 'diff.mids')[c(1, 1:bins, bins)],
	 c(0,attr(x, 'diff.dist'),0), col=col,
	 border=do.call(rgb, as.list(c(col2rgb(col)/256, 0.5))));
   }
}

enve.TRIBS.merge <- function
   ### Merges two `enve.TRIBS` objects generated from the same objects at
   ### different subsampling levels.
   (x,
   ### First `enve.TRIBS` object.
   y
   ### Second `enve.TRIBS` object.
   ){
   # Check consistency
   if(attr(x,'distance') != attr(y,'distance'))
      stop('Total distances in objects are different.');
   if(any(attr(x,'points') != attr(y,'points')))
      stop('Points in objects are different.');
   if(attr(x,'spaceSize') != attr(y,'spaceSize'))
      stop('Space size in objects are different.');
   if(attr(x,'selSize') != attr(y,'selSize'))
      stop('Selection size in objects are different.');
   if(attr(x,'dimensions') != attr(y,'dimensions'))
      stop('Dimensions in objects are different.');
   if(nrow(attr(x,'distances')) != nrow(attr(y,'distances')))
      stop('Replicates in objects are different.');
   # Merge
   a <- attr(x,'subsamples');
   b <- attr(y,'subsamples');
   o <- order(c(a,b));
   o <- o[!duplicated(c(a,b)[o])] ;
   d <- cbind(attr(x,'distances'), attr(y,'distances'))[, o] ;
   z <- new('enve.TRIBS',
      distance=attr(x,'distance'), points=attr(x,'points'),
      distances=d, spaceSize=attr(x,'spaceSize'),
      selSize=attr(x,'selSize'), dimensions=attr(x,'dimensions'),
      subsamples=c(a,b)[o], call=match.call());
   return(z) ;
   ### Returns an `enve.TRIBS` object.
}

#==============> Define core functions
enve.tribs.test <- function
   ### Estimates the empirical difference between all the distances in a set of
   ### objects and a subset, together with its statistical significance.
   (dist,
   ### Distances as `dist` object.
   selection,
   ### Selection defining the subset.
   bins=50,
   ### Number of bins to evaluate in the range of distances.
   ...
   ### Any other parameters supported by `enve.tribs`, except `subsamples`.
   ){
   s.tribs <- enve.tribs(dist, selection, subsamples=c(0,1), ...);
   a.tribs <- enve.tribs(dist,
      subsamples=c(0,attr(s.tribs, 'selSize')/attr(s.tribs, 'spaceSize')), ...);
   s.dist <- attr(s.tribs, 'distances')[, 2];
   a.dist <- attr(a.tribs, 'distances')[, 2];
   range <- range(c(s.dist, a.dist));
   a.f <- hist(a.dist, breaks=seq(range[1], range[2], length.out=bins),
      plot=FALSE);
   s.f <- hist(s.dist, breaks=seq(range[1], range[2], length.out=bins),
      plot=FALSE);
   zp.f <- c(); zz.f <- 0; zn.f <- c();
   p.x <- a.f$counts/sum(a.f$counts);
   p.y <- s.f$counts/sum(s.f$counts);
   for(z in 1:length(a.f$mids)){
      zn.f[z] <- 0;
      zz.f <- 0;
      zp.f[z] <- 0;
      for(k in 1:length(a.f$mids)){
         if(z < k){
	    zp.f[z] <- zp.f[z] + p.x[k]*p.y[k-z];
	    zn.f[z] <- zn.f[z] + p.x[k-z]*p.y[k];
	 }
	 zz.f <- zz.f + p.x[k]*p.y[k];
      }
   }
   return(new('enve.TRIBStest',
      pval.gt=sum(c(zz.f, zp.f)), pval.lt=sum(c(zz.f, zn.f)),
      all.dist=p.x, sel.dist=p.y, diff.dist=c(rev(zn.f), zz.f, zp.f),
      dist.mids=a.f$mids,
      diff.mids=seq(diff(range(a.f$mids)), -diff(range(a.f$mids)),
      length.out=1+2*length(a.f$mids)),
      call=match.call()));
   ### Returns an `enve.TRIBStest` object.
}

enve.tribs <- function
   ### Subsample any objects in "distance space" to reduce the effect of
   ### sample-clustering.  This function was originally designed to subsample
   ### genomes in "phylogenetic distance space", a clear case of strong
   ### clustering bias in sampling, by Luis M. Rodriguez-R and Michael R
   ### Weigand.
   (dist,
   ### Distances as a `dist` object.
   selection=labels(dist),
   ### Objects to include in the subsample. By default, all objects are
   ### selected.
   replicates=1000,
   ### Number of replications per point
   summary.fx=median,
   ### Function to summarize the distance distributions in a given replicate. By
   ### default, the median distance is estimated.
   dist.method='euclidean',
   ### Distance method between random points and samples in the transformed
   ### space. See `dist`.
   subsamples=seq(0,1,by=0.01),
   ### Subsampling fractions
   dimensions=ceiling(length(selection)*0.05),
   ### Dimensions to use in the NMDS. By default, 5% of the selection length.
   metaMDS.opts=list(),
   ### Any additional options to pass to metaMDS, as `list`.
   threads=2,
   ### Number of threads to use.
   verbosity=1,
   ### Verbosity. Use 0 to run quietly, increase for additional information.
   points,
   ### Optional. If passed, the MDS step is skipped and this object is used
   ### instead.  It can be the `$points` slot of class `metaMDS` (from `vegan`).
   ### It must be a matrix or matrix-coercible object, with samples as rows and
   ### dimensions as columns.
   pre.tribs
   ### Optional. If passed, the points are recovered from this object (except if
   ### `points` is also passed. This should be an `enve.TRIBS` object estimated
   ### on the same objects (the selection is unimportant).
   ){
   if(!is(dist, 'dist'))
      stop('`dist` parameter must be a `dist` object.');
   # 1. NMDS
   if(missing(points)){
      if(missing(pre.tribs)){
	 if(verbosity > 0)
	    cat('===[ Estimating NMDS ]\n');
	 if(!suppressPackageStartupMessages(
	    requireNamespace("vegan", quietly=TRUE)))
	       stop('Unavailable required package: `vegan`.');
	 mds.args <- c(metaMDS.opts, list(comm=dist, k=dimensions,
	    trace=verbosity));
	 points <- do.call(vegan::metaMDS, mds.args)$points;
      }else{
	 points <- attr(pre.tribs, 'points');
	 dimensions <- ncol(points);
      }
   }else{
      points <- as.matrix(points);
      dimensions <- ncol(points);
   }
   # 2. Pad ranges
   if(verbosity > 0) cat('===[ Padding ranges ]\n');
   dots <- matrix(NA, nrow=nrow(points), ncol=dimensions,
      dimnames=list(rownames(points), 1:dimensions));
   selection <- selection[!is.na(match(selection, rownames(dots)))];
   for(dim in 1:dimensions){
      dimRange <- range(points[,dim]) +
	 c(-1,1)*diff(range(points[,1]))/length(selection);
      dots[, dim] <- (points[,dim]-dimRange[1])/diff(dimRange);
   }
   # 3. Select points and summarize distances
   if(verbosity > 0) cat('===[ Sub-sampling ]\n');
   distances <- matrix(NA, nrow=replicates, ncol=length(subsamples),
      dimnames=list(1:replicates, as.character(subsamples)));
   cl <- makeCluster(threads);
   for(frx in subsamples){
      if(verbosity > 1) cat('Sub-sampling at ',(frx*100),'%\n',sep='');
      distances[, as.character(frx)] = parSapply(cl, 1:replicates, enve.__tribs,
	 frx, match(selection, rownames(dots)), dimensions, dots, dist.method,
	 summary.fx, dist);
   }
   stopCluster(cl);
   # 4. Build object and return
   return(new('enve.TRIBS',
      distance=do.call(summary.fx, list(as.matrix(dist)[selection, selection])),
      points=points, distances=distances, spaceSize=nrow(points),
      selSize=length(selection), dimensions=dimensions, subsamples=subsamples,
      call=match.call()));
   ### Returns an `enve.TRIBS` object.
}

enve.__tribs <- function
   ### Internal ancilliary function (see `enve.tribs`).
   (rep, frx, selection, dimensions, dots, dist.method, summary.fx, dist){
   sample <- c();
   if(frx==0) return(0);
   for(point in 1:round(frx*length(selection))){
      rand.point <- runif(dimensions);
      closest.dot <- '';
      closest.dist <- Inf;
      for(dot in selection){
	 dot.dist <- as.numeric(dist(matrix(c(rand.point, dots[dot,]), nrow=2,
	    byrow=TRUE), method=dist.method));
	 if(dot.dist < closest.dist){
	    closest.dot <- dot;
	    closest.dist <- dot.dist;
	 }
      }
      sample <- c(sample, closest.dot);
   }
   return( do.call(summary.fx, list(as.matrix(dist)[sample, sample])) );
}


