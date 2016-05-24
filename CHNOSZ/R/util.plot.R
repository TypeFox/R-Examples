# CHNOSZ/util.plot.R
# Functions to create and modify plots

thermo.plot.new <- function(xlim,ylim,xlab,ylab,cex=par('cex'),mar=NULL,lwd=par('lwd'),side=c(1,2,3,4),
  mgp=c(1.5,0.3,0),cex.axis=par('cex'),col=par('col'),yline=NULL,axs='i',do.box=TRUE,ticks=NULL,
  las=1,xline=NULL) {
  # start a new plot with some customized settings
  thermo <- get("thermo")
  # 20120523 store the old par in thermo$opar
  if(is.null(thermo$opar)) {
    thermo$opar <- par(no.readonly=TRUE)
    assign("thermo", thermo, "CHNOSZ")
  }
  # 20091108 changed argument name from 'ticks' to 'side' but
  # keep 'ticks' for backward compatibility
  if(!is.null(ticks)) side <- ticks 
  # 20090324 mar handling: NULL - a default setting; NA - par's setting
  # 20090413 changed mar of top side from 2 to 2.5
  if(is.null(mar)) mar <- c(3,3.5,2.5,1) else if(is.na(mar[1])) mar <- par('mar')
  par(mar=mar,mgp=mgp,tcl=0.3,las=las,xaxs=axs,yaxs=axs,cex=cex,lwd=lwd,col=col,fg=col)
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  if(do.box) box()
  # labels
  if(is.null(xline)) xline <- mgp[1]
  thermo.axis(xlab,side=1,line=xline,cex=cex.axis,lwd=NULL)
  if(is.null(yline)) yline <- mgp[1]
  thermo.axis(ylab,side=2,line=yline,cex=cex.axis,lwd=NULL)
  # (optional) tick marks
  if(1 %in% side) thermo.axis(NULL,side=1,lwd=lwd)
  if(2 %in% side) thermo.axis(NULL,side=2,lwd=lwd)
  if(3 %in% side) thermo.axis(NULL,side=3,lwd=lwd)
  if(4 %in% side) thermo.axis(NULL,side=4,lwd=lwd)
}

thermo.axis <- function(lab='x-axis',side=1,line=1.5,cex=par('cex'),lwd=par('lwd'),T=NULL,col=par('col')) {
  # if T isn't NULL, looks like we want make a second
  # oxidation scale corresponding to one already plotted.
  # e.g.,  Eh-pe, Eh-logfO2, or logfO2-Eh
  if(!is.null(T)) {
    usr <- par('usr')
    if(side %in% c(1,3)) lim <- usr[1:2] else lim <- usr[3:4]
    if(length(grep('pe',lab)) > 0) {
      lim <- convert(lim,'pe',T=T)
    } else if(length(grep('O2',lab)) > 0) {
      lim <- convert(lim,'logfO2',T=T)
    } else if(length(grep('Eh',lab)) > 0) {
      lim <- convert(lim,'E0',T=T)
    }
    if(side %in% c(1,3)) usr[1:2] <- lim else usr[3:4] <- lim
    opar <- par(usr=usr)
  }
  if(!is.null(lwd)) {
    ## plot major tick marks and numeric labels
    do.label <- TRUE
    if(missing(cex) & side %in% c(3,4) & is.null(T)) do.label <- FALSE
    at <- axis(side,labels=do.label,tick=TRUE,lwd=lwd,col=col,col.axis=col) 
    ## plot minor tick marks
    # the distance between major tick marks
    da <- abs(diff(at[1:2]))
    # distance between minor tick marks
    di <- da / 4
    if(da %% 2 | !(da %% 10)) di <- da / 5
    # number of minor tick marks
    if(side %in% c(1,3)) {
      ii <- c(1,2) 
      myasp <- par('xaxp')
    } else {
      ii <- c(3,4)
      myasp <- par('yaxp')
    }
    myusr <- par('usr')[ii]
    daxis <- abs(diff(myusr))
    nt <- daxis / di + 1
    ## if nt isn't an integer, it probably
    ## means the axis limits don't correspond
    ## to major tick marks (expect problems)
    ##at <- seq(myusr[1],myusr[2],length.out=nt)
    # start from (bottom/left) of axis?
    bl <- 1
    #if(myasp[2]==myusr[2]) bl <- 2
    # is forward direction (top/right)?
    tr <- 1
    if(xor(myusr[2] < myusr[1] , bl==2)) tr <- -1
    #at <- myusr[bl] + tr * di * seq(0:(nt-1))
    # well all of that doesn't work in a lot of cases,
    # where none of the axis limits correspond to
    # major tick marks. perhaps the following will work
    at <- myusr[1] + tr * di * (0:(nt-1))
    # apply an offset
    axt <- axTicks(side)[1]
    daxt <- (axt - myusr[1])/di
    daxt <- (daxt-round(daxt))*di
    at <- at + daxt
    tcl <- par('tcl') * 0.5
    axis(side,labels=FALSE,tick=TRUE,lwd=lwd,col=col,col.axis=col,at=at,tcl=tcl)
  }

  # rotate labels on side axes
  if(side %in% c(2,4)) las <- 0 else las <- 1
  if(!is.null(lab)) mtext(lab,side=side,line=line,cex=cex,las=las)
  # reset limits if we were plotting a second axis
  if(!is.null(T)) par(opar)
}

label.plot <- function(x, xfrac=0.05, yfrac=0.95, paren=FALSE, italic=FALSE, ...) {
  # make a text label e.g., "(a)" in the corner of a plot
  # xfrac, yfrac: fraction of axis where to put label (default top right)
  # paren: put a parenthesis around the text, and italicize it?
  if(italic) x <- substitute(italic(a), list(a=x))
  if(paren) x <- substitute(group('(',a,')'), list(a=x))
  if(italic | paren) x <- as.expression(x)
  pu <- par('usr')
  text(pu[1]+xfrac*(pu[2]-pu[1]), pu[3]+yfrac*(pu[4]-pu[3]), labels=x, ...)
}

usrfig <- function() {
  # function to get the figure limits in user coordinates
  # get plot limits in user coordinates (usr) and as fraction [0,1] of figure region (plt)
  xusr <- par('usr')[1:2]; yusr <- par('usr')[3:4]
  xplt <- par('plt')[1:2]; yplt <- par('plt')[3:4]
  # linear model to calculate figure limits in user coordinates
  xlm <- lm(xusr ~ xplt); ylm <- lm(yusr ~ yplt)
  xfig <- predict.lm(xlm, data.frame(xplt=c(0, 1)))
  yfig <- predict.lm(ylm, data.frame(yplt=c(0, 1)))
  return(list(x=xfig, y=yfig))
}

label.figure <- function(x, xfrac=0.05, yfrac=0.95, paren=FALSE, italic=FALSE, ...) {
  # function to add labels outside of the plot region  20151020
  f <- usrfig()
  # similar to label.plot(), except we have to set xpd=TRUE here
  opar <- par(xpd=NA)
  if(italic) x <- substitute(italic(a), list(a=x))
  if(paren) x <- substitute(group('(',a,')'), list(a=x))
  if(italic | paren) x <- as.expression(x)
  text(f$x[1]+xfrac*(f$x[2]-f$x[1]), f$y[1]+yfrac*(f$y[2]-f$y[1]), labels=x, ...)
  par(opar)
}

water.lines <- function(xaxis='pH', yaxis='Eh', T=298.15, P='Psat', which=c('oxidation','reduction'),
  logaH2O=0, lty=2, lwd=1, col=par('fg'), xpoints=NULL, O2state="gas") {
  # draw water stability limits
  # if we're on an Eh-pH diagram, or logfO2-pH diagram,
  # or logfO2-T or Eh-T
  # calculate them exactly (nicer looking lines), otherwise 
  # (TODO) add them using affinity() and diagram()
  
  # get the x and y limits from the current plot
  pu <- par('usr')
  xlim <- pu[1:2]
  ylim <- pu[3:4]
  # exact lines
  # warning: Eh calculations are reliable only at a single T
  if(xaxis=="O2" | (xaxis=='pH' & (yaxis=='Eh' | yaxis=='O2' | yaxis=="pe"))) {
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c("H2O", "O2", "H2"), c(-1, 0.5, 1), c("liq", O2state, "gas"), T=T, P=P, convert=FALSE)$out$logK 
      # this is logfO2 if O2state=="gas", or logaO2 if O2state=="aq"
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(xaxis=='O2') abline(v=logfO2,lty=lty,lwd=lwd,col=col) 
      else if(yaxis=='O2') abline(h=logfO2,lty=lty,lwd=lwd,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,lwd=lwd,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,lwd=lwd,col=col)
    }
    if('oxidation' %in% which) {
      logfO2 <- 0
      logK <- subcrt(c("O2", "O2"), c(-1, 1), c("gas", O2state), T=T, P=P, convert=FALSE)$out$logK 
      # this is logfO2 if O2state=="gas", or logaO2 if O2state=="aq"
      logfO2 <- logfO2 + logK
      if(xaxis=='O2') abline(v=logfO2,lty=lty,lwd=lwd,col=col) 
      if(yaxis=='O2') abline(h=logfO2,lty=lty,lwd=lwd,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,lwd=lwd,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,lwd=lwd,col=col)
    }
  } else if(xaxis %in% c('T','P') & yaxis %in% c('Eh','O2') ) {
    #if(xaxis=='T') if(is.null(xpoints)) xpoints <- T
    # 20090212 get T values from plot limits
    # TODO: make this work for T on y-axis too
    if(xaxis=='T') {
      if(missing(T)) {
        xpoints <- seq(xlim[1],xlim[2],length.out=100)
        T <- envert(xpoints,"K")
      }
    }
    if(xaxis=='P') if(is.null(xpoints)) xpoints <- P
    if('oxidation' %in% which) {
      logfO2 <- rep(0,length(xpoints))
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,lwd=lwd,col=col)
      else lines(xpoints,logfO2,lty=lty,lwd=lwd,col=col)
    }
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c('H2O','oxygen','hydrogen'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,lwd=lwd,col=col)
      else lines(xpoints,logfO2,lty=lty,lwd=lwd,col=col)
    }
  } else {
    # inexact lines
    #
  }
}

mtitle <- function(main, line=0, ...) {
  # make a possibly multi-line plot title 
  # useful for including expressions on multiple lines 
  # 'line' is the margin line of the last (bottom) line of the title
  l <- length(main)
  for(i in 1:l) mtext(main[i], line=line+l-i, ...)
}


residualsplot <- function(residuals,property="Cp",model="big") {
  # a barchart useful for plotting residuals 
  # (see 'xadditivy' vignette)
  # plot horizontal bars stacked with largest on top
  # put them in order of increasing size
  residuals <- residuals[order(abs(residuals))]
  bp <- barplot(residuals,horiz=TRUE,axisnames=FALSE)
  # to plot names, first find the sign of the residual
  ispos <- residuals > 0
  # put some padding around the names
  nr <- paste("",names(residuals),"")
  # align left or right depending on sign of residual
  text(0,bp[ispos],nr[ispos],adj=1)
  text(0,bp[!ispos],nr[!ispos],adj=0)
  # assemble stuff for the title
  title(main=paste("residuals in",property,"using",model,"model"))
  # add axis label
  mtext(axis.label(property),1,3)
}
