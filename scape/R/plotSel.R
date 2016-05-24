plotSel <- function(model, together=FALSE, series=NULL, sex=NULL, axes=TRUE, legend="bottom", main="", xlab="", ylab="",
                    cex.main=1.2, cex.legend=1, cex.lab=1, cex.axis=0.8, cex.strip=0.8, col.strip="gray95",
                    strip=strip.custom(bg=col.strip), las=1, tck=0, tick.number=5, lty.grid=3, col.grid="gray", pch="m",
                    cex.points=1, col.points="black", lty.lines=1, lwd.lines=4, col.lines=c("red","blue"), plot=TRUE,
                    ...)
{
  ## 1  Define functions
  panel.each <- function(x, y, subscripts, maturity, col.lines.vector, ...)
  {
    panel.grid(h=-1, v=-1, lty=lty.grid, col=col.grid)
    panel.points(maturity$Age, maturity$P, col=col.points, ...)
    panel.lines(x, y, col=col.lines.vector[subscripts], ...)
  }
  panel.together <- function(x, y, subscripts, maturity, ...)
  {
    panel.grid(h=-1, v=-1, lty=lty.grid, col=col.grid)
    panel.points(maturity$Age, maturity$P, col=col.points, ...)
    panel.superpose(x, y, type="l", subscripts=subscripts, col=col.lines, ...)
  }

  ## 2  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))

  ## 3  Prepare data (extract, rearrange, filter)
  x <- model$Sel
  if(is.null(series))
    series <- unique(x$Series)
  if(is.null(sex))
    sex <- unique(x$Sex)
  ok.series <- x$Series %in% series; if(!any(ok.series)) stop("Please check if the 'series' argument is correct")
  ok.sex <- x$Sex %in% sex; if(!any(ok.sex)) stop("Please check if the 'sex' argument is correct")
  x <- x[ok.series & ok.sex,]
  if(is.numeric(x$Series))
    x$Series <- factor(paste("Series", x$Series))
  mat <- x[x$Series=="Maturity",]
  sel <- x[x$Series!="Maturity",]
  sel$Series <- factor(as.character(sel$Series))  # update levels

  ## 4  Prepare plot (set pars, vectorize args, create list args)
  nseries <- length(unique(sel$Series))
  lty.lines <- rep(lty.lines, length.out=max(2,nseries))
  lwd.lines <- rep(lwd.lines, length.out=max(2,nseries))
  col.lines <- rep(col.lines, length.out=max(2,nseries))
  mymain <- list(label=main, cex=cex.main)
  myxlab <- list(label=xlab, cex=cex.lab)
  myylab <- list(label=ylab, cex=cex.lab)
  myrot <- switch(as.character(las), "0"=list(x=list(rot=0),y=list(rot=90)), "1"=list(x=list(rot=0),y=list(rot=0)),
                  "2"=list(x=list(rot=90),y=list(rot=0)), "3"=list(x=list(rot=90),y=list(rot=90)))
  myscales <- c(list(draw=axes,cex=cex.axis,tck=tck,tick.number=tick.number), myrot)
  mystrip <- strip.custom(bg=col.strip)
  mytext <- list(cex=cex.strip)
  mykey <- list(space=legend, text=list(lab=levels(sel$Series),cex=cex.legend),
                lines=list(lty=lty.lines,lwd=lwd.lines,col=col.lines))

  ## 5  Create trellis object
  if(!together)
  {
    graph <- xyplot(P~Age|Series*Sex, data=sel, panel=panel.each, maturity=mat, as.table=TRUE,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    pch=pch, cex=cex.points, col.points=col.points, lty=lty.lines, lwd=lwd.lines,
                    col.lines.vector=col.lines[factor(x$Sex)], ...)
  }
  else
  {
    graph <- xyplot(P~Age|Sex, data=sel, groups=sel$Series, panel=panel.together, maturity=mat,
                    main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip, par.strip.text=mytext,
                    key=mykey,
                    pch=pch, cex=cex.points, col.points=col.points, lty=lty.lines, lwd=lwd.lines, col.line=col.lines,
                    ...)
  }
  if(is.list(graph$x.limits))                                                            # set xlim=0,max(ages)&ylim=0,1
  {
    graph$x.limits <- lapply(graph$x.limits, function(x){x<-c(0,max(x$Age));return(x)})  # multi-panel plot
    graph$y.limits <- lapply(graph$y.limits, function(y){y<-c(-0.005,1.005);return(y)})  # multi-panel plot
  }
  else
  {
    graph$x.limits <- c(0, max(x$Age))                                                   # single-panel plot
    graph$y.limits <- c(-0.005,1.005)                                                    # single-panel plot
  }
  ## 6  Finish
  if(plot)
  {
    print(graph)
    invisible(x)
  }
  else
  {
    invisible(graph)
  }
}
