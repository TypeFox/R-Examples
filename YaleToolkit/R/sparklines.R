"sparklines" <-
function(ss,
                             times = NULL,
			     overlap = FALSE,
                             yscale = NULL,
                             buffer = unit(0, 'lines'),
                             buffer.pars = NULL,
                             IQR = NULL,
                             ptopts = NULL,
                             yaxis = TRUE,
                             xaxis = "exterior",
                             labeled.points = NULL,
                             point.labels = NULL, label.just = c(1.2, 0.5),
                             frame.pars = NULL,
                             line.pars = gpar(lwd = 1),
                             outer.margin = unit(c(5,4,4,2), 'lines'),
                             outer.margin.pars = NULL,
                             main = NULL,
                             sub = NULL,
                             xlab = NULL, #need to implement
                             ylab = NULL,
                             lcol=NULL,
                             new = TRUE) {


  if (is.null(ss) && !is.null(outer.margin.pars)) {
    grid.rect(gp = outer.margin.pars)
    return()
  }

  if(!is.data.frame(ss))
    stop('ss is not a data frame')
  if(is.null(times))
    times <- 1:nrow(ss)
  if(length(times) != nrow(ss))
    stop('the length of times is not the same as the number of observations')
  if (is.list(yscale)) {
    if (length(yscale)!=length(ss)) stop("incorrect dimensions of ss and yscale")
    yscales <- yscale
  } else {
    yscales <- vector('list', length = length(ss))
    if (is.vector(yscale) && length(yscale) == 2)
      for (j in 1:length(ss)) yscales[[j]] <- yscale
    else {
      if (is.null(yscale)) {
        yscale <- vector('list', length = length(ss)) 
        for (j in 1:length(ss)) yscales[[j]] <- range(ss[,j], na.rm = TRUE)
      }
    }
  }
  if (is.null(lcol)) lcol <- rep(1, length(ss))

  if (new) grid.newpage()

  if(!is.null(ss) && !is.null(outer.margin.pars))
    grid.rect(gp = outer.margin.pars)

  # A buffer if desired (default is not)...
  subvp <- viewport(x = outer.margin[2], y = outer.margin[1],
                    width = unit(1, 'npc') - outer.margin[2] - outer.margin[4],
            height = unit(1, 'npc') - outer.margin[1] - outer.margin[3],
            just = c('left', 'bottom'))
  pushViewport(subvp)

  if (!is.null(main)) grid.text(main, x=unit(0.5, "npc"), 
                                y=unit(1, "npc")+unit(1.5, "lines"),
                                gp=gpar(fontface=2))
  if(!overlap){
    panel.layout <- viewport(layout = grid.layout(length(ss),1),
                             xscale=range(times, na.rm = TRUE))
    pushViewport(panel.layout)
  
    if (!is.null(xaxis) & xaxis=="exterior") {
      grid.xaxis()
      if (!is.null(xlab)){
        grid.text(label = xlab, x = unit(0.5, 'npc'),
                  y = unit(-3, "lines"))
      }
      xaxis <- FALSE
      xlab <- NULL
    }
    ## Loop through sparklines, printing them not overlapped
    for(i in 1:length(ss)){
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = i,
                          yscale=yscales[[i]]))
      sparkline(s=ss[,i], times=times, ylim=yscales[[i]],
                    buffer = buffer, new=FALSE,
                    ptopts = ptopts, #c(1,3,8),
                    frame.pars = frame.pars,
                    buffer.pars = buffer.pars,
                    yaxis = yaxis, ylab=ylab[i], 
                    xaxis = xaxis, xlab=xlab[i],
                    sub = sub[i],
                    IQR = IQR, line.pars=gpar(col=lcol[i]))
      ## The following 2 (or 3) pops are relocated from sparkline()
      popViewport(1) # pops subvp
      if (xaxis==TRUE || xaxis=='exterior') popViewport(1) # pops outervp
  	  popViewport(1) # pops outer.margins
      ## end of relocation from sparkline()
      popViewport() # pop this cell in layout
    }
    popViewport(1) # pop the panel.layout
  }else{ # if overlap is TRUE

    if (!is.null(xaxis) & xaxis=="exterior") {
      grid.xaxis()
      if (!is.null(xlab)){
        grid.text(label = xlab, x = unit(0.5, 'npc'),
                  y = unit(-3, "lines"))
      }
      xaxis <- FALSE
      xlab <- NULL
    }

    ## Loop through sparklines, printing them overlapped
    ## THERE'S A PROBLEM if xaxis is not specified...why?

    if(length(unique(yscales)) > 1){
      warning("y-scales are not the same; are you really sure you want to plot all sparklines on the same y-axis?")
    }

    for(i in 1:length(ss)){
      pushViewport(viewport(y = 0, height = 1,
                   yscale = yscales[[i]], just = "bottom",
                   default.units = "npc"))
      sparkline(s=ss[,i], times=times, ylim=yscales[[i]],
                buffer = buffer, new=FALSE,
                ptopts = ptopts,
                frame.pars = frame.pars,
                buffer.pars = buffer.pars,
                yaxis = yaxis, ylab=ylab[i], 
                xaxis = xaxis, xlab=xlab[i],
                sub = sub[i],
                IQR = NULL,
                line.pars=gpar(col=lcol[i]))
	}
  }
  
  popViewport(1) # pop subvp

} # End of function sparklines

