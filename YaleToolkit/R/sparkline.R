"sparkline" <-
function(s,
                            times = NULL,
                            ylim = NULL,
                            buffer = unit(0, 'lines'),
                            margins = NULL,
                            IQR = NULL,
                            yaxis = FALSE,
                            xaxis = FALSE,
                            ptopts = list(points = NULL,
                                          labels = NULL,
                                          labels.ch = NULL,
                                          gp = NULL,
                                          just = NULL,
                                          pch = NULL),
                            margin.pars = NULL,
                            buffer.pars = NULL,
                            frame.pars = NULL,
                            line.pars = gpar(lwd = 1),
                            main = NULL,
                            sub = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            new = TRUE) {

  ###################################################################
  # Sanity checking, option processing, and then create the viewport:

  if(is.ts(s)){
    times <- seq(attributes(s)$tsp[1], attributes(s)$tsp[2],
                 by = attributes(s)$tsp[3])
  }
  
  if(!is.null(times) && (length(times)!=length(s) || any(is.na(times))))
    warning("inconsistency in times; ignoring times.")

  if(is.null(times) || length(times)!=length(s) || any(is.na(times)))
    times <- time(s)
  xlim <- range(times)

  if (is.vector(ptopts) && !is.list(ptopts)) {
    if (is.numeric(ptopts))
      ptopts <- list(points = ptopts, labels = ptopts,
                     labels.ch = as.character(ptopts))
    else
      ptopts <- list(points = NULL, labels = ptopts, labels.ch = NULL)
  }
  if (is.null(ptopts$gp)) ptopts$gp <- gpar(col = 'red')
  if (is.null(ptopts$pch)) ptopts$pch <- 19
  if (is.null(ptopts$just)) ptopts$just <- c(-0.3, 0.5)
  if (any(ptopts$labels=='first.last')) {
    ptopts$labels <- as.numeric(ptopts$labels[ptopts$labels!='first.last'])
    ptopts$labels <- c(ptopts$labels, 1, length(s))
    ptopts$points <- c(ptopts$points, 1, length(s))
    ptopts$labels.ch <- c(ptopts$labels.ch, format(s[c(1, length(s))], digits = 2))
  }
  if (any(ptopts$labels=='min.max')) {
    ptopts$labels <- as.numeric(ptopts$labels[ptopts$labels!='min.max'])
    ptopts$labels <- c(ptopts$labels,
                       (1:length(s))[s==min(s, na.rm=TRUE) |
                                     s==max(s, na.rm=TRUE)])
    ptopts$points <- c(ptopts$points,
                       (1:length(s))[s==min(s, na.rm=TRUE) |
                                     s==max(s, na.rm=TRUE)])
    tmplabels <- rep(NA, length(ptopts$labels))
    tmplabels[s[ptopts$labels]==min(s, na.rm=TRUE)] <-
                                 format(min(s, na.rm=TRUE), digits = 2)
    tmplabels[s[ptopts$labels]==max(s, na.rm=TRUE)] <-
                                 format(max(s, na.rm=TRUE), digits = 2)    
    if (is.null(ptopts$labels.ch)) ptopts$labels.ch <- tmplabels
    else ptopts$labels.ch <- c(ptopts$labels.ch, tmplabels)
  }

  if (is.null(ylim)) ylim <- range(s, na.rm = TRUE)
  else if (length(ylim) != 2) stop('ylim is the wrong length')  
  
  if (new){
    grid.newpage()
    if(is.null(margins)) margins <- unit(c(0.02,0.02,0.02,0.02), 'npc')
    if(!is.null(sub)){   # CHECK THIS!!!
      margins <- unit.c(unit(c(5,4,4), 'lines'), unit(1, 'strwidth', sub) + unit(3, "lines"))
    }
  } else {
    if(is.null(margins)) margins <- unit(c(0,0,0,0), 'lines')
  }

  if (!is.null(margin.pars)) grid.rect(gp = margin.pars)
  outer.margins <- viewport(x = margins[2], y = margins[1],
                    width = unit(1, 'npc') - margins[2] - margins[4],
            height = unit(1, 'npc') - margins[1] - margins[3],
            just = c('left', 'bottom'))
  pushViewport(outer.margins)

  if (!is.null(main)) grid.text(main, x=unit(0.5, "npc"), y=unit(1, "npc") + unit(1.5, "lines"),
                                gp=gpar(fontface='bold'))
  if (!is.null(sub)) grid.text(sub, y=unit(0.5, "npc"), x=unit(1, "npc") + unit(1.5, "lines"),
                                just = 'left', gp=gpar(fontface='plain'))

  if (xaxis==TRUE || xaxis=='exterior') { 
    outervp = viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       width = unit(1, "npc"), height = unit(1, "npc"),
                       xscale = xlim, yscale = ylim, just = "center")
    pushViewport(outervp)
    grid.xaxis()
  }
  grid.text(xlab, unit(0.5, "npc"), unit(-3, "lines"))
  if (!is.null(buffer.pars)) grid.rect(gp = buffer.pars)

  subvp <- viewport(x = unit(0.5, 'npc'),
                    y = unit(0.5, 'npc'),
                    width = unit(1, 'npc'),
                    height = unit(1, 'npc') - 2 * buffer,
                    xscale = unit(xlim, "native"),
                    yscale = unit(ylim, "native"),
                    just = 'center')
  pushViewport(subvp)
  if (!is.null(frame.pars)) grid.rect(gp = frame.pars)
  if (yaxis) {
    grid.yaxis()
    grid.text(ylab, unit(-3.5, "lines"), unit(0.5, "npc"), rot=90)
  }
  if (xaxis=='interior') grid.xaxis()

  ##########################################################
  # Now that we have the viewport, do the interesting stuff:

  if (!is.null(IQR)) {
    b <- boxplot(s, plot = FALSE)$stats
    grid.rect(x = unit(0.5, 'npc'),
      width = unit(1, 'npc'),
      y = unit(b[2,1], 'native'),
      height = unit(abs(b[4,1] - b[2,1]), 'native'), gp = IQR,
              just = c("center", "bottom"))
  }
  
  # This is where the time series is actually plotted
  grid.lines(times, s, default.units = 'native', gp = line.pars)
  
  # Plotting points:
  if (!is.null(ptopts$points)) {
    grid.points(times[ptopts$points], s[ptopts$points],
           gp = ptopts$gp, default.units = 'native', pch = ptopts$pch)
  }

  # Labelling points:
  if (!is.null(ptopts$labels)) {
    if (is.null(ptopts$labels.ch)) ptopts$labels.ch = as.character(ptopts$labels)
    grid.text(label = ptopts$labels.ch,
              x = times[ptopts$labels],
              y = s[ptopts$labels],
              default.units = 'native', just = ptopts$just)
  }
  
  ## These are commented out because I think we want the fct to return
  ##  with the actual plotting vp active, but maybe this will break
  ##  something....moved to sparklines()
  #popViewport(1) # pops subvp
  #if (xaxis==TRUE || xaxis=='exterior') popViewport(1) # pops outervp
  #popViewport(1) # pops outer.margins

} # End of function sparkline
