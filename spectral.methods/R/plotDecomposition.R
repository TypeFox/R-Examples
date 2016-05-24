plotDecomposition = function(
  ##title<< Plot the results of a SSA decomposition
  decomp.results    ##<< list: object to plot: results from a call to filterTSeriesSSA
   ,time.vector=c() ##<< R Date/time object: optional time vector used to label x axes.
   , n.magnif=3     ##<< amount of positions where to plot the different magnifications.
                    ##   This results in the amount of columns in the plot
   )
  ##description<< This function visualizes the results from a call to filterTSeriesSSA
  ##value<<
  ## Nothing is returned.
{
  if (!(class(decomp.results)=='list' &
        (sum(names(decomp.results)==c('dec.series','borders'))==2)))
    stop('decomp.results does not seem to be the result from a call to filterTSeriesSSA!')

  dec.series      <-  decomp.results$dec.series
  borders         <-  decomp.results$borders
  data.decomposed <-  decomp.results$data.decomposed
  if (length(time.vector)==0)
    time.vector=1:dim(dec.series)[2]
  n.bands      <- dim(borders)[1]
  n.dpts       <- dim(dec.series)[2]
  n.cycle.plot <- c(3,12)

  x.range.band     <- integer(length=n.bands)
  mean.period.band <- apply(borders,MARGIN=1,mean)
  period.band.t    <- n.dpts
  x.range.band[mean.period.band==Inf] <- n.dpts
  step=1
  while(is.element(0,x.range.band)& step < n.bands+1) {
    bands.scale.t <- mean.period.band < (period.band.t*n.cycle.plot[1]) &
      mean.period.band > (period.band.t/n.cycle.plot[2])
    x.range.band[x.range.band==0 & bands.scale.t] <- period.band.t
    if (!sum(x.range.band==0)==0)
      period.band.t <- max(mean.period.band[x.range.band==0])
    step=step+1
  }
  n.rows.plot <-length(unique(x.range.band[!x.range.band==n.dpts]))
  n.cols.plot <- n.magnif
  x.range.rows<- sort(unique(x.range.band[!x.range.band==n.dpts]),decreasing=TRUE)

  ## do plots
  layout(matrix(c(rep(1,times=n.cols.plot),2:(n.cols.plot*n.rows.plot+1)),
                ncol=n.cols.plot,byrow=TRUE))
  par(mar=c(2,0,0,0),oma=c(2,4,2,4),tcl=0.2,mgp=c(0.5,0.5,0))
  ind.first    <- sort(sample(1:(dim(data.decomposed$dec.series)[2]-
                                 max(x.range.rows*mean(n.cycle.plot))),n.cols.plot))
  series.recstr<- colSums(dec.series)
  yrange.orig <- diff(range(series.recstr) )
  plot(time.vector,series.recstr,pch='.',cex=2,xlab='',ylab='')
  legend('topright',merge=TRUE,lty=c(NA,rep(1,n.bands)),pch=c('.',rep(NA,times=n.bands)),
         col=1:(n.bands+1),legend=c('rec. series',apply(borders,MARGIN=1,function(x)paste(x,collapse='-'))),
         cex=2,ncol=2,pt.cex=4)
  bands.origscale <- which(x.range.band==n.dpts)
  offset.scale.points=(min(series.recstr)-(min(dec.series[bands.origscale,])))
  for (h in bands.origscale)
    points(time.vector,dec.series[h,]+offset.scale.points,
           type='l',col=h+1,lwd=2)
  abline(v=c(ind.first),lty=2)
  if (length(bands.origscale)>0) {
    ylim.points <- par()$usr[3:4]-offset.scale.points
    par(new=TRUE)
    plot(1,1,ylim= ylim.points,type='n',xaxt='n',xlab='',ylab='',yaxt='n')
    axis(4)
    }
  
  for(i in 1:n.rows.plot) { 
    for(j in 1:n.cols.plot) {
      ind.plot      <- seq(from=ind.first[j],length.out=floor(x.range.rows[i]*mean(n.cycle.plot)))
      cex=8-floor(log(length(ind.plot),10))
      series.plot.t <- which(x.range.band==x.range.rows[i])
      plot(time.vector[ind.plot],series.recstr[ind.plot],ylim=range(series.recstr),yaxt='n',pch='.',cex=cex,xlab='')
      par(new=TRUE)
      plot(range(ind.plot),mean(dec.series[series.plot.t,])+(c(-0.5,0.5)*yrange.orig),type='n',yaxt='n',xlab='')
      if (j==1) {
        axis(2)
      } else if (j == n.cols.plot) {
        axis(4)
      }        
      for(k in series.plot.t)
        points(time.vector[ind.plot],dec.series[k,ind.plot],type='l',xaxt='n',yaxt='n',col=k+1)
    }
  }
  mtext(outer=TRUE,side=1,line=0,text='time index')
  mtext(outer=TRUE,side=2,line=2,text='scale reconstructed series')
  mtext(outer=TRUE,side=4,line=2,text='scale decomposed parts')
}
