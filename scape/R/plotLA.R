plotLA <- function(model, together=FALSE, sex=NULL, axes=TRUE, same.limits=TRUE, between=list(x=axes,y=axes), ylim=NULL,
                   bands=1, main="", xlab="", ylab="", cex.main=1.2, cex.lab=1, cex.axis=0.8, cex.strip=0.8,
                   col.strip="gray95", strip=strip.custom(bg=col.strip), las=1, tck=0, tick.number=5, lty.grid=3,
                   col.grid="gray", pch=16, cex.points=0.5, col.points="black", lty.lines=1, lwd.lines=4,
                   col.lines=c("red","blue"), lty.bands=2*!together, lwd.bands=1, col.bands="black", plot=TRUE, ...)
{
  ## 1  Define functions
  panel.each <- function(x, y, subscripts, col.points, col.lines, col.bands, ...)  # obs, fit, and bands in sex panels
  {
    panel.grid(h=-1, v=-1, lty=lty.grid, col=col.grid)
    panel.superpose.2(x, y, subscripts=subscripts, col.symbol=col.points[subscripts],
                      col.line=c(col.lines[subscripts[1]],col.bands[subscripts]), ...)
  }
  panel.together <- function(x, y, z, subscripts, pch, cex.points, col.points, col.lines, ...)  # obs, fit, and bands
  {
    panel.grid(h=-1, v=-1, lty=lty.grid, col=col.grid)
    one <- z==levels(z)[1]
    two <- z==levels(z)[2]
    panel.superpose.2(x[one], y[one], subscripts=subscripts[one], pch=pch[1], cex=cex.points[1],
                      col.symbol=col.points[1], lty=c(lty.lines[1],lty.bands[1]), lwd=c(lwd.lines[1],lwd.bands[1]),
                      col.line=col.lines[1], ...)
    panel.superpose.2(x[two], y[two], subscripts=subscripts[two], pch=pch[2], cex=cex.points[2],
                      col.symbol=col.points[2], lty=c(lty.lines[2],lty.bands[1]), lwd=c(lwd.lines[2],lwd.bands[1]),
                      col.line=col.lines[2], ...)
  }
  ## 2  Parse args
  if(class(model) != "scape")
    stop("The 'model' argument should be a scape object, not ", class(model))
  relation <- if(same.limits) "same" else "free"

  ## 3  Prepare data (extract, rearrange, filter)
  x <- model$LA
  obs <- x[!is.na(x$Obs),]
  fit <- aggregate(list(Fit=x$Fit,CV=x$CV), list(Sex=x$Sex,Age=x$Age), mean)
  fit$Age <- as.integer(as.character(fit$Age))
  fit <- fit[order(fit$Sex,fit$Age),]
  x <- data.frame(Sex=c(as.character(obs$Sex),as.character(rep(fit$Sex,3))),
                  Age=c(obs$Age,rep(fit$Age,3)),
                  ObsFit=c(rep("Obs",nrow(obs)),rep("Fit",nrow(fit)),rep("Upper",nrow(fit)),rep("Lower",nrow(fit))),
                  Length=c(obs$Obs,fit$Fit,fit$Fit*exp(bands*fit$CV),fit$Fit/exp(bands*fit$CV)))
  if(is.null(sex))
    sex <- unique(x$Sex)
  ok.sex <- x$Sex %in% sex; if(!any(ok.sex)) stop("Please check if the 'sex' argument is correct")
  x <- x[ok.sex,]
  if(!bands)
    x <- x[x$ObsFit!="Upper"|x$ObsFit!="Lower",]
  nsexes <- length(unique(x$Sex))

  ## 4  Prepare plot (set pars, vectorize args, create list args)
  pch <- rep(pch, length.out=2)
  cex.points <- rep(cex.points, length.out=2)
  col.points <- rep(col.points, length.out=2)
  lty.lines <- rep(lty.lines, length.out=2)
  lwd.lines <- rep(lwd.lines, length.out=2)
  col.lines <- rep(col.lines, length.out=2)
  col.bands <- rep(col.bands, length.out=2)
  mymain <- list(label=main, cex=cex.main)
  myxlab <- list(label=xlab, cex=cex.lab)
  myylab <- list(label=ylab, cex=cex.lab)
  myrot <- switch(as.character(las), "0"=list(x=list(rot=0),y=list(rot=90)), "1"=list(x=list(rot=0),y=list(rot=0)),
                  "2"=list(x=list(rot=90),y=list(rot=0)), "3"=list(x=list(rot=90),y=list(rot=90)))
  myscales <- c(list(draw=axes,relation=relation,cex=cex.axis,tck=tck,tick.number=tick.number), myrot)
  mystrip <- strip.custom(bg=col.strip)
  mytext <- list(cex=cex.strip)

  ## 5  Create trellis object
  if(is.null(ylim))
    ylim <- c(0, 1.04*max(x$Length))
  if(nsexes==2 && together)
  {
    graph <- xyplot(Length~Age, data=x, groups=x$ObsFit, z=factor(x$Sex), panel=panel.together, type=c("l","l","p"),
                    ylim=ylim, main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip,
                    par.strip.text=mytext,
                    pch=pch, cex=cex.points, col.points=col.points, col.lines=col.lines, ...)
  }
  else
  {
    graph <- xyplot(Length~Age|Sex, data=x, groups=x$ObsFit, panel=panel.each, type=c("l","l","p"), as.table=TRUE,
                    between=between,
                    ylim=ylim, main=mymain, xlab=myxlab, ylab=myylab, scales=myscales, strip=strip,
                    par.strip.text=mytext,
                    pch=pch, cex=cex.points, col.points=col.points[factor(x$Sex)], lty=c(lty.lines[1],lty.bands[1]),
                    lwd=c(lwd.lines[1],lwd.bands[1]), col.lines=col.lines[factor(x$Sex)],
                    col.bands=col.bands[factor(x$Sex)], ...)
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
