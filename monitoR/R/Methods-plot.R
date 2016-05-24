# For plotting templates
# Modified: 6 Sept 2015

setMethod('plot', signature(x='TemplateList', y='ANY'), 
  function(
    x,                             # Complete template list
    which.one=names(x@templates),  # If template list is provided, which template(s) should be used?
    click=FALSE,                   # Set to FALSE for no interaction
    ask=if(length(which.one)>1) TRUE else FALSE, # Set to FALSE for no pause between plots
    spec.col=gray.3(),             # Color palette for spectrogram
    on.col='#FFA50075',            # Color for on points in binary templates
    off.col='#0000FF75',           # Color for off points in binary templates
    pt.col='#FFA50075',            # Color for correlation template points
    line.col='black'               #
  ) {
  
    # Set mai and get oldask
    oldpar <- par(mai=c(1.02, 0.82, 0.82, 0.82))
    oldask <- par(ask=par('ask'))
    on.exit(par(c(oldpar, oldask)))
    
    # Template loop
    for(i in which.one) {
      template <- x@templates[[i]]
  
      # Pull out clip name
      clip.path <- template@clip.path
      file.ext <- tolower(gsub(".*\\.", "", clip.path))
      if(file.ext == 'wav') clip <- tuneR::readWave(filename=clip.path) else if(file.ext == 'mp3') clip.path <- readMP3(filename=clip) else stop('File extension must be wav or mp3, got ', file.ext)
  
      # Get time and frequency limits
      # t.lim is for cutting the clip, and needs to include samples (i.e., wave@left) that are hanging at the end when ovlp > 0
      # The next line actually should return one extra point
      t.lim <- template@first.t.bin + c(0, template@duration + template@t.step*100/(100-template@ovlp))
      frq.lim <- template@frq.lim + c(-1, 1)*template@frq.step/2
  
      wave <- cutWave(clip, from=t.lim[1], to=t.lim[2]) 
  
      # Fourier transform
      fspec <- spectro(wave=wave, wl=template@wl, ovlp=template@ovlp, wn=template@wn)
      # Filter amplitudes 
      t.bins <- fspec$time
      n.t.bins <- length(t.bins)
      which.t.bins <- 1:n.t.bins
      which.frq.bins <- which(fspec$freq >= frq.lim[1] & fspec$freq <= frq.lim[2])
      frq.bins <- fspec$freq[which.frq.bins]
      n.frq.bins <- length(frq.bins)
      amp <- fspec$amp[which.frq.bins, ]
      t.step <- t.bins[2]-t.bins[1]
      frq.step <- frq.bins[2]-frq.bins[1]
  
      # Create plot
      image(x=which.t.bins, y=which.frq.bins, t(amp), col=spec.col, xlab='Time (s)', ylab='Frequency (kHz)', las=1, useRaster=TRUE, axes=FALSE, las=1)
      t.bin.ticks <- pretty(t.bins+t.lim[1], n=6);axis(1, at=(t.bin.ticks-t.lim[1])/t.step, labels=t.bin.ticks)
      frq.bin.ticks <- pretty(frq.bins, n=6);axis(2, at=frq.bin.ticks/frq.step, labels=frq.bin.ticks, las=1)
      axis(3);axis(4, las=1);box()
      mtext(paste('Template', i), 3,line=3, cex=1.2)
   
      # Plot template points
      if(class(x) == 'binTemplateList') {
  
        pt.on <- template@pt.on
        pt.off <- template@pt.off
  
        pt.on[, 'frq'] <- pt.on[, 'frq'] - min(which.frq.bins) + 1
        pt.off[, 'frq'] <- pt.off[, 'frq'] - min(which.frq.bins) + 1
        bin.amp <- 0*amp
        bin.amp[pt.on[, c(2, 1)]] <- 1
        bin.amp[pt.off[, c(2, 1)]] <- 2
  
        image(x=which.t.bins, y=which.frq.bins, t(bin.amp), zlim=c(0, 2), col=c('transparent', on.col, off.col), add=TRUE)
  
      } else if (class(x) == 'corTemplateList') {
  
        pts <- template@pts
        pts[, 'frq'] <- pts[, 'frq'] - min(which.frq.bins) + 1
        bin.amp <- 0*amp
        bin.amp[pts[, c(2, 1)]] <- 1
  
        image(x=which.t.bins, y=which.frq.bins, t(bin.amp), zlim=c(0, 1), col=c('transparent', pt.col), add=TRUE)
  
      } else stop('Template list class not recognized.')
  
      # Start loop for identifying plot locations 
      if(click) {
        i <- 0
        click.pts <- NULL
        message('Identify plot locations with left mouse click. To exit, right click.')
        pos <- NULL
        while(!is.null(pos)|i == 0) {
          i <- i + 1
          pos <- locator(n=1)
          if(!is.null(pos)) {
            pos <- lapply(pos, round)
            abline(h=pos$y, col=line.col)
            abline(v=pos$x, col=line.col)
            text(pos$x, pos$y, i,col='red')
  
            # locator returned bin positions, so time and frequency values need to be determined from indexing
            t.pos <- t.lim[1] + t.bins[pos$x] + t.step # + t.step there because first time bin is zero. . .IS THIS A PROBLEM ELSEWHERE?
            frq.pos <- frq.bins[pos$y - min(which.frq.bins) + 1]
            amp.pos <- amp[pos$y - min(which.frq.bins) + 1, pos$x]
            click.pts <- rbind(click.pts, c(t.pos, frq.pos, amp.pos))
  
            x.mid <- mean(par()$usr[1:2])
            y.lim <- par()$usr[3:4]
  
            text(x.mid, y.lim[2]-i/20*diff(y.lim), paste(i, '.time=', signif(t.pos, 3), ',frq=', signif(frq.pos, 3), ',amp=', signif(amp.pos, 3), sep=''), cex=0.8, col=line.col)
          }
        } 
      }
      par(ask=ask)
    }

    if(click) {
      colnames(click.pts) <- c('t', 'frq', 'amp')
      invisible(click.pts)
    }

  }
)


# templateScores

setMethod('plot', signature(x='detectionList', y='ANY'), 
  function(
  x,  
  flim=c(0, 12),                 # Frequency limits for the spectrogram.
  scorelim,                     # Plot limits for scores
  which.one=names(x@templates), # Name(s) of templates to plot
  box=TRUE,                     # Set to FALSE to surpress boxes in spectrogram showing hits
  spec.col=gray.2(),            # Color palette for spectrogram
  t.each=30,                    # Time shown for each individual plot (s)
  hit.marker='lines',           # Markers for hits in score plot
  color=c('red', 'blue', 'green', 'orange', 'purple', 'pink', 'darkgreen', 'turquoise', 'royalblue', 'orchid4', 'brown', 'salmon2'), # Colors for individual templates
  legend=TRUE,                  # Set to FALSE to surpress legend
  all.peaks=FALSE,              # Set to TRUE to indicate locations of all peaks
  ask=if(dev.list() == 2) TRUE else FALSE 
  ) {

    survey <- x@survey
    t.survey <- length(survey@left)/survey@samp.rate 
    n.plots <- ceiling(t.survey/t.each)
    t.start <- 1:n.plots*t.each - t.each
    if(n.plots == 1) t.each <- t.survey
    t.end <- t.start + t.each
    t.end[t.end>t.survey] <- t.survey
    t.start[n.plots] <- t.end[n.plots] - t.each # Adjust start of last plot back so it has the same length as the others

    # Pull out spectrogram data from scores object
    # Based on first template
    amp <- x@survey.data[[1]]$amp
    t.bins <- x@survey.data[[1]]$t.bins
    frq.bins <- x@survey.data[[1]]$frq.bins
 
    # Sort out colors for lines and boxes
    names.t <- names(x@templates)
    n.templates <- length(names.t)
    color <- c(rep(color, n.templates %/% length(color)), color[1:n.templates%%length(color)])
    names(color) <- names.t

    # Get scorelim
    if(missing(scorelim)) {
      upr <- 0
      for(i in seq(length(x@scores))) {
        upr <- max(upr, x@scores[[i]]$score)
      }
      scorelim <- c(0, upr)
    }

    oldpar <- par(mar=c(1, 4,1, 1), oma=c(6, 0,0, 0), mfrow=c(2, 1))
    oldask <- par(ask=par('ask'))
    on.exit(par(c(oldpar, oldask)))

    # Loop through time windows, plotting a spectrogram for each time
    for(i in 1:length(t.start)) {

      message(paste(t.start[i], 'to', t.end[i], 'seconds'))
    
      times <- t.bins[t.bins>=t.start[i] & t.bins<=t.end[i]]
      amp.clip <- amp[, t.bins %in% times]
      image(x=times, y=frq.bins, t(amp.clip), ylim=flim, col=spec.col, xlab='', ylab='Frequency (kHz)', xaxt='n', las=1)
 
      # Loop through templates and add boxes around detections
      for(j in which.one) {
        template <- x@templates[[j]]
        if(all.peaks) pks <- x@peaks[[j]] else pks <- x@detections[[j]]
        pks.clip <- pks[pks$time + template@duration >= t.start[i] & pks$time - template@duration <= t.end[i], ]

        if(box & nrow(pks.clip)>0) {
          for(k in 1:nrow(pks.clip)) { 
            xleft <- pks.clip$time[k] - template@duration/2
            xright <- pks.clip$time[k] + template@duration/2
            ylwr <- template@frq.lim[1]
            yupr <- template@frq.lim[2]
            polygon(x=c(xleft, xleft, xright, xright), y=c(ylwr, yupr, yupr, ylwr), border=color[j], lwd=1)
          }
        }
      }

      # Make plot of scores. Can't sort out xlab for some reason.
      plot(NULL, xlim=c(t.start[i], t.end[i]), ylim=scorelim, xlab='', ylab='Score', type='n', xaxs='i', las=1, mgp=c(3, 1,0))
      mtext("Time (s or min:sec)", 1,2.5, outer=TRUE)

      # Add x axis as mm:ss
      xaxp.sec <- par('xaxp')
      labs.sec <- seq(xaxp.sec[1], xaxp.sec[2], length.out=xaxp.sec[3]+1)
      labs.mmss <- paste(sprintf('%02d', labs.sec%/%60), ':', sprintf('%02d', labs.sec%%60), sep='')
      axis(1, at=labs.sec, labels=labs.mmss, mgp=c(3, 1.9, 0))

      if(legend) legend('topright', which.one, lty=1, col=color[which.one], cex=0.7)

      # Loop through templates 
      for(j in which.one) {
        template <- x@templates[[j]]
        score <- x@scores[[j]]         # score output from sccDetect. The correlation coefficients within which hits were found.
        if(all.peaks) pks <- x@peaks[[j]] else pks <- x@detections[[j]]
        cutoff <- template@score.cutoff      # If given, will plot a horizontal line at the correlation coefficient cutoff.

        score.clip <- score[score$time>=t.start[i] & score$time<=t.end[i], ]
        pks.clip <- pks[pks$time + template@duration >= t.start[i] & pks$time - template@duration <= t.end[i], ]

        lines(score.clip$time, score.clip$score, col=color[j])
        if(hit.marker == 'points') points(pks.clip$time, pks.clip$score, col=color[j]) else
          if(hit.marker == 'lines') abline(v=pks.clip$time, col=color[j]) 
        if(is.vector(cutoff)) abline(h=cutoff, lty=2, col=color[j])
      }
      par(ask=ask)
    }
  
  }
)  


