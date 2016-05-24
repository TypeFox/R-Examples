.lc.main <- 
function(y, type,
       col.line, col.area, col.stroke, col.fill, shape.pts,
       col.grid, col.box, col.bg, cex.axis, col.axis,
       rotate.values, offset, xy.ticks,
       line.width, xlab, ylab, main, sub, cex,
       time.start, time.by, time.reverse, 
       center.line, show.runs, quiet, ...) {


  # scale for regular R or RStudio
  adj <- .RSadj(bubble.size=NULL, cex.axis)
  size.axis <- adj$size.axis
  size.lab <- adj$size.lab

  if (!is.numeric(y)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The variable must be numeric.\n\n")
  }

  if (!is.null(type)) if (type != "p" && type != "l" && type != "b") { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Option 'type' can only be \"p\" for points,\n",
        "  \"l\" for line or \"b\" for both.\n\n")
  }

  if (!is.null(time.start) && is.null(time.by)) {
     cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Specified  time.start  so also need  time.by.\n\n")
  }

  nrows <- length(y)
  pt.size <- ifelse (is.null(cex), 0.8, cex)

  # get variable label and axis labels if they exist
  gl <- .getlabels(ylab=ylab, main=main, cex.lab=getOption("lab.size"))
  y.name <- gl$yn;  y.lbl <- gl$yl;  y.lab <- gl$yb
  main.lab <- gl$mb
  sub.lab <- gl$sb
  cex.lab <- gl$cex.lab

  # count and remove missing data
  n <- sum(!is.na(y))
  n.miss <- sum(is.na(y))
  y <- na.omit(y)

  x <- numeric(length = 0)
  if (time.reverse) {  # use x as a temporary place holder
    for (i in 1:nrows) x[i] <- y[nrows+1-i]
    y <- x
  }

  if (is.null(time.start) && class(y) != "ts") {
    x.lab <- ifelse (is.null(xlab), "index", xlab)
    x <- 1:length(y)  # ordinal position of each value on x-axis
  }
  else   {  # time.start date specified or a ts
    x.lab <- ifelse (is.null(xlab), "", xlab)      
    if (class(y) == "ts") {
      time.start <- paste(start(y)[1], "/", start(y)[2], "/01", sep="") 
      frq <- frequency(y)
      if (frq == 1) time.by <- "year"
      if (frq == 4) time.by <- "3 months"
      if (frq == 12) time.by <- "month"
      if (frq == 52) time.by <- "week"
      if (frq == 365) time.by <- "year"
    }
    date.seq <- seq.Date(as.Date(time.start), by=time.by, length.out=nrows)
    x <- date.seq  # dates on x-axis
  }

  if (center.line == "default") {  # by default display center.line only if many runs
    m <- mean(y, na.rm=TRUE)
    n.change <- 0
    for (i in 1:(length(y)-1)) if ((y[i+1]>m) != (y[i]>m)) n.change <- n.change+1 
    if (n.change/(length(y)-1) < .15) center.line <- "off" else center.line <- "median"
  }

  # fill ts chart 
  #if (!is.null(time.start) && is.null(col.area))
    #col.area <- getOption("col.fill.bar")

  if (is.null(type)) 
    if (is.null(col.area) || col.area == "transparent") type <- "b" 
    else type <- "l"


  # set point size       
  if (type == "b" && is.null(cex))  
    if (nrows < 50) pt.size <- 0.8 else pt.size <- .9 - 0.002*nrows
  if (pt.size < 0) pt.size <- 0
  
  digits.d <- .max.dd(y) + 1
  options(digits.d=digits.d)

  if (is.null(main)) {
    orig.params <- par(no.readonly=TRUE)
    on.exit(par(orig.params))
    par(mar=c(4,4,2,2)+0.1)
  }

  # plot setup
  plot(x, y, type="n", axes=FALSE, ann=FALSE, ...)

  if (xy.ticks){
    if (is.null(time.start) && class(x)!="ts") 
     .axes(x.lvl=NULL, y.lvl=NULL, axTicks(1), axTicks(2),
        par("usr")[1], par("usr")[3], size.axis, col.axis,
        rotate.values, offset, ...)
    else {
      axis.Date(1, x, cex.axis=size.axis, col.axis=col.axis, ...)
      #lbl.dt <- as.Date(axTicks(1), origin = "1970-01-01")
      #axis.Date(1, x, labels=FALSE, tck=-.01, ...)
      #text(x=lbl.dt, y=par("usr")[3], labels=lbl.dt,
           #pos=1, xpd=TRUE, cex=cex.axis, col=col.axis)
      axis(2, at=axTicks(2), labels=FALSE, tck=-.01, ...)
      dec.d <- .getdigits(round(axTicks(2),6),1) - 1
      text(x=par("usr")[1], y=axTicks(2), labels=.fmt(axTicks(2),dec.d),
           pos=2, xpd=TRUE, cex=size.axis, col=col.axis)
    }
  }

  # axis labels
  max.lbl <- max(nchar(axTicks(2)))
  .axlabs(x.lab, y.lab, main.lab, sub.lab, max.lbl, 
          xy.ticks=TRUE, offset=offset, cex.lab=size.lab, ...) 

  # colored plotting area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col=col.bg, border=col.box)

  # grid lines
  vx <- pretty(c(usr[1],usr[2]))
  vy <- pretty(c(usr[3],usr[4]))
  abline(v=seq(vx[1],vx[length(vx)],vx[2]-vx[1]), col=col.grid, lwd=.5)
  abline(h=seq(vy[1],vy[length(vy)],vy[2]-vy[1]), col=col.grid, lwd=.5)

  # fill area under curve
  if (!is.null(col.area)  && !is.null(col.stroke)) {
    if (col.area=="transparent"  &&  col.stroke=="transparent")
      col.area <- NULL
    else if (type != "p"  &&  is.null(col.stroke))
      col.stroke <- col.line
  }
    if (!is.null(col.area)) 
      polygon(c(x[1],x,x[length(x)]), c(min(y),y,min(y)),
              col=col.area, border=col.stroke)

  # the plot
  if (type == "l" || type == "b") {
    lines(as.numeric(x),y, col=col.line, lwd=line.width, ...)
  }
  if (type == "p" || type == "b") {
    points(x,y, col=col.stroke, pch=shape.pts, bg=col.fill, cex=pt.size, ...)
  }

  # plot center line
  if (center.line != "off") {
    if (center.line == "mean") {
      m.y <- mean(y)
      lbl <- "mean"
      lbl.cat <- "mean:"
    }
    else if (center.line == "median") {
      m.y <- median(y)
      lbl <- "medn"
      lbl.cat <- "median:"
    }
    else if (center.line == "zero") {
      m.y <- 0
      lbl <- ""
      lbl.cat <- "median:"
    }

    abline(h=m.y, col="gray50", lty="dashed")
    mtext(lbl, side=4, cex=.9, col="gray50", las=2, at=m.y, line=0.1)
    if (center.line == "zero") m.y <- median(y) 

    gl <- .getlabels()
    x.name <- gl$xn; x.lbl <- gl$xl;
    y.name <- gl$yn; y.lbl <- gl$yl
    if (!quiet) {
      ttlns <- .title2(x.name, y.name, x.lbl, y.lbl, TRUE)
      class(ttlns) <- "out_piece"
      output <- list(out_title=ttlns)
      class(output) <- "out_all"
      print(output)
    }


    # analyze runs
    if (!quiet) {
      cat("\n")
      cat("n:", n, "\n")
      cat("missing:", n.miss, "\n")
      cat(lbl.cat, round(m.y,digits.d), "\n")
      cat("\n")
      .dash(12); cat("Run Analysis\n"); .dash(12)
      run <- integer(length=0)  # length of ith run in run[i]
      n.runs <- 1  # total number of runs
      run[n.runs] <- 1
      line.out <- "    1"
      cat("\n")
      for (i in 2:length(y)) {
        if (y[i] != m.y) {  # throw out values that equal m.y
          if (sign(y[i]-m.y) != sign(y[i-1]-m.y)) {  # new run
            if (show.runs) {
              if (n.runs < 10) buf <- "  " else buf <- " "
              cat("size=", run[n.runs], "  Run", buf, n.runs, ":",
                  line.out, "\n", sep="")
            }
            line.out <- ""
            n.runs <- n.runs + 1
            run[n.runs] <- 0
          }
          run[n.runs] <- run[n.runs] + 1
          buf <- ifelse (i < 10, "  ", " ")
          line.out <- paste(line.out, buf, i)
        }
      }
      cat("size=", run[n.runs], "  Run", buf, n.runs, ":", line.out,
          "\n", sep="")
      eq.ctr <- which(y==m.y)
      cat("\nTotal number of runs:", n.runs, "\n")
      txt <- "Total number of values that do not equal the "
      cat(txt, lbl.cat, " ", length(y)-length(eq.ctr), "\n", sep="")
      if (length(eq.ctr) != 0) {
        if (show.runs) {
          cat("\nValues ignored that equal the", lbl.cat, "\n")
          for (i in 1:length(eq.ctr))
            cat("    #", eq.ctr[i], " ", y[eq.ctr[i]], sep="", "\n")
          cat("Total number of values ignored:", length(eq.ctr), "\n")
        }
      }
      else 
        cat("Total number of values ignored that equal the", lbl.cat, 
            length(eq.ctr), "\n")
    }
  }

  if (!quiet) cat("\n")

}

