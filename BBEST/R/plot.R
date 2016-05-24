################################################################
# test.signal(x, lambda, sigma, x.delta, knots.n, peaks.widthRange,  peaks.n)
#
# returns:      Generates a random function with smooth background, returning 
#               separate curves for the background, the peaks, and the noise.
#   list() with the following elements:
#     $x:       x-values for this function (same as this function received)
#     $curves:  full-res curves (data.frame() with following columns):
#       $y:       total curve (bkg + peaks + noise)
#       $bkg:     background contribution
#       $SB:    zero
#     $noise:  IID Gaussian noise
#     $knots:  spline knots (data.frame() with following columns):
#       $x:      x-values of knots
#       $y:      y-values of knots
# arguments:
#   x:                 numeric vector of x-values
#   lambda:            maximum half-height of any peak
#   sigma:             noise level
#   x.delta:           minimum spacing between consecutive spline knots
#   knots.n:           number of spline knots for background (less means smoother)
#   peaks.widthRange:  range in peak widths
#   peaks.n:           number of peaks to add
test.signal <- function(x, lambda, sigma, x.delta, knots.n, peaks.widthRange,  peaks.n) {
  # 1.  Generate peaks.
  # Assume mean height roughly half peak height (hence "2/lambda"):
  peaksAmplitude <- rexp(n=peaks.n, rate=(2 / lambda))
  # random peaks width 
  peaksWidth <- runif(n=peaks.n, min=peaks.widthRange[1], max=peaks.widthRange[2])
  # separate peak position by 2 sigma
  peaksCenter <- runif(n=peaks.n,
        min=min(x) + (2 * peaksWidth),
        max=max(x) - (2 * peaksWidth))
  # Assume Gaussian peak shape
  N <- length(x)
  peaks <- 0
  for(i in 1:N){
    peaks[i] <- 0
	for(j in 1:peaks.n)
	  peaks[i] <- peaks[i] + peaksAmplitude[j]*exp(-0.5 * ((x[i] - peaksCenter[j]) / peaksWidth[j]) ^ 2)
  } 
  
  # 2. Add noise
  noise <- rnorm(sd=sigma, n=N)  # add Poisson noise?

  # 3. Construct background 
  knots.x <- 0
  for(i in 1:knots.n){
    knots.x[i] <-  runif(n=1, min=min(x), max=max(x)-(knots.n-1)*x.delta)	
  }
  knots.x <- sort(knots.x)
  knots.x <- knots.x +  ((1:knots.n) - 1) * x.delta # this procedure was checked and works correctly
  knots.y <- rnorm(n=knots.n, sd=0.3)
  bkg <- get.bkg(x=x, knots.x=knots.x, knots.y=knots.y)
  if(min(bkg) < 0){
	knots.y <- knots.y + abs(min(bkg))
    bkg <- bkg + abs(min(bkg))
  }
  # 4. put everything together
  noisyPeaks <- peaks + noise
  signal <- noisyPeaks + bkg
    
  return (list(x=x, y=signal, sigma=rep(sigma, length(x)), 
               SB=rep(0, length(x)), lambda=rep(lambda, length(x)),
			   knots=list(x=knots.x, y=knots.y), bkg=bkg))
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# (Cookbook for R)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#########################################################
#mPlot.results(fit.results, label.x, label.y){
#
# returns:       plots background estimate and corrected signal
# arguments   
#   fit.results: return value of do.fit
#   label.x,  
#   label.y:     plot labels
mPlot.results <- function(fit.results, label.x="x", label.y="y", xlim=NA, ylim=NA){
 # require(ggplot2)
	if(length(fit.results$uncrt) == 1 && is.na(fit.results$pars)){
	  dat <- list(x=fit.results$x, y=fit.results$curves$y, SB=fit.results$curves$SB, 
              lambda=fit.results$fit.details$lambda, sigma=fit.results$fit.details$sigma)
	  if(!is.null(fit.results$curves$corr))
	    dat$y <- dat$y + fit.results$curves$corr
	  fit.results$uncrt <- get.hess(dat, fit.results$knots$x, fit.results$knots$y, Gr=NA, r=seq(0, 2, 0.01), p.bkg=.5)
	}
  
  if(!is.na(fit.results$pars)){
    fit.results$uncrt <- list()  
    fit.results$uncrt$stdev <- rep(0, length(fit.results$x))
  }
  
  if(any(is.na(ylim)))
    ylim <- c(min(fit.results$curves$y), max(fit.results$curves$y))    
  if(any(is.na(xlim)))
    xlim <- c(min(fit.results$x), max(fit.results$x))

  if(any(is.na(fit.results$knots$x)) || any(is.na(fit.results$knots$y))){
      fit.results$knots$x <- 0
      fit.results$knots$y <- 0  
    }
  #to avoid NOTEs in package check:
  x=y=bkg=stdev=signal=variable=value=NULL
	curves <- data.frame(x = fit.results$x, y = fit.results$curves$y-fit.results$curves$SB, 
	                     bkg = fit.results$curves$bkg)				 
	melted.curves <- reshape2::melt(curves, id.vars="x")
	ribbon <- data.frame(x = fit.results$x, bkg = fit.results$curves$bkg, stdev = fit.results$uncrt$stdev,
	                     signal = fit.results$curves$y-fit.results$curves$bkg)
	knots <- data.frame(x=fit.results$knots$x, y=fit.results$knots$y)
	p1<- (h <- ggplot2::ggplot(melted.curves, aes(x=x, y=value))
     + geom_line(aes(colour=variable, group=variable)) 
     + scale_color_manual(name="Plot legend", values = c("black", "red", "brown"), labels = c("experimental data", "background guess (+/-2sd)", "uncertainty interval"))
     + geom_ribbon(data=ribbon, aes(ymin=bkg-2*stdev, ymax=bkg+2*stdev, y=bkg), fill = "brown", alpha=0.5)
     + geom_point(data=knots, aes(x=x, y=y), colour="red", size=2)
     + xlab(label.x) 
     + ylab(label.y) 
     + xlim(xlim)
     + ylim(ylim)
     + ggtitle("Background Estimation")
	)
#  print(p1)

	curves2 <- data.frame(x = fit.results$x, signal = fit.results$curves$y-fit.results$curves$bkg, SB=fit.results$curves$SB)
	melted.curves2 <- reshape2::melt(curves2, id.vars="x")
  ylim[1] <- ylim[1] + min(fit.results$curves$SB-fit.results$curves$bkg)
  ylim[2] <- ylim[2] + max(fit.results$curves$SB-fit.results$curves$bkg)
	p2<- (h2 <- ggplot2::ggplot(melted.curves2, aes(x=x, y=value))
     + geom_line(aes(colour=variable, group=variable)) 
     + scale_color_manual(name="Plot legend", values = c("blue", "green", "gray20"), labels = c("corrected signal", "coherent baseline"))
     + geom_ribbon(data=ribbon, aes(ymin=signal-2*stdev, ymax=signal+2*stdev, y=signal), fill = "gray20", alpha=0.5)
     + xlab(label.x) 
     + ylab(label.y) 
     + xlim(xlim)
#     + ylim(ylim)     
     + ggtitle("Corrected Signal")
	)
    	

	multiplot(p1, p2, cols=1)
}

#########################################################
# calc.Gr(fit.results, rho.0, r.min, r.max, Q.min, Q.max, nsd)
#
# returns:       the PDF function. Also plots PDF and corresponding confidence interval
# arguments   
#   fit.results: return value of do.fit
#   rho.0:       atomic number density for the material
#   r.min,       plot G(r) form r.min to r.max
#   r.max:  
#   Q.min,       cut S(Q) at the interval between Q.min 
#   Q.max:       and Q.max
#   nsd:         number of standard devitions to plot the uncertainty


calc.Gr <- function(fit.results, rho.0, plot=TRUE, r.min=0, r.max=5, dr=0.01, Q.min=NA, Q.max=NA, nsd=2, gr.compare=NA){
  r <- seq(r.min, r.max, dr)
  SQ <- fit.results$curves$y-fit.results$curves$bkg
  Q <- fit.results$x
  if(is.na(Q.max)) Q.max <- max(Q)
  if(is.na(Q.min)) Q.min <- Q.max*.95
  ind.max <- which(abs(Q-Q.max)==min(abs(Q-Q.max)))
  ind.min <- which(abs(Q-Q.min)==min(abs(Q-Q.min)))
  cut <- which(abs(SQ[ind.min:ind.max]-1) ==min(abs(SQ[ind.min:ind.max]-1) ))[1] + ind.min - 1
   
  dat <- list(x=Q[1:cut], y=fit.results$curves$y[1:cut], SB=fit.results$curves$SB[1:cut], 
              lambda=fit.results$fit.details$lambda[1:cut], sigma=fit.results$fit.details$sigma[1:cut])
  if(!is.na(fit.results$fit.details$Gr[1])){
    cat("Recalculating G(r) information... \n")
    dat <- set.Gr(dat, r1=r, rho.0=rho.0, type1="gaussianNoise")			 
  }
  else
    dat$Gr <- NA
  

  cat("Calculating standard deviation... \n")
  if(!is.null(fit.results$curves$corr))
	  dat$SB <- dat$SB - fit.results$curves$corr[1:cut]
  
  if(is.na(fit.results$pars)){  
    knots.x <- fit.results$knots$x
    knots.y <- fit.results$knots$y
    stdev <- get.hess(dat, knots.x, knots.y, Gr=dat$Gr, r=r, p.bkg=.5)$stdev.r
  }
  else
    stdev <- rep(0, length(r))
    
  if(!is.null(fit.results$curves$corr))
	  dat$SB <- dat$SB + fit.results$curves$corr[1:cut]
  
  cat("Calculating Pair Distribution Function... \n")
  gr <- sineFT(f.Q=dat$y-fit.results$curves$bkg[1:cut]-1, Q=dat$x, r=r)

  stdev <- stdev*nsd
  if(plot==TRUE)
    fplot.Gr(r=r, gr=gr, stdev=stdev, rho.0=rho.0, nsd=nsd, gr.compare=gr.compare)
    
  return(list(r=r, gr=gr, stdev=stdev/nsd))
}

fplot.Gr <- function(r, gr, stdev, rho.0, nsd=2, gr.compare=NA, xlim=NA, ylim=NA, title="corrected G(r)"){
 # require(ggplot2)
  cat("Plotting... \n")
  if(any(is.na(ylim)))
    ylim=c(min(gr)*1.1-abs(max(stdev)), max(gr)*1.1+abs(max(stdev)))
  else{
    ylim[1] <- ylim[1]*1.1 - abs(max(stdev))
    ylim[2] <- ylim[2]*1.1 + abs(max(stdev))
  }
  if(any(is.na(xlim)))
    xlim=c(min(r), max(r))
  #to avoid NOTEs in package check:
  x=y=variable=value=NULL
  if(!is.na(gr.compare[1])){
    curves <- data.frame(x = r, y = gr, gr.compare=gr.compare, l = -4*pi*rho.0*r)	 
    vals = c("black", "red", "darkblue", "blue")
    labs = c(title, "G(r) to compare", paste("uncertainty interval (+/-",nsd,"sd)", sep=""))
  }
  else{  
    curves <- data.frame(x = r, y = gr, l = -4*pi*rho.0*r)	 
    vals = c("black", "darkblue", "blue")
    labs = c(title, paste("uncertainty interval (+/-",nsd,"sd)", sep=""))
  }
  melted.curves <- reshape2::melt(curves, id.vars="x")
  ribbon <- data.frame(x = r, y = gr, stdev = stdev)
  options(warn=-1)
  p1<-(h <- ggplot2::ggplot(melted.curves, aes(x=x, y=value))
    + geom_line(aes(colour=variable, group=variable))
    + geom_ribbon(data=ribbon, aes(ymin=y-stdev, ymax=y+stdev, y=y), fill = "blue", alpha=0.5)
    + scale_color_manual(name="Plot legend", values = vals, labels = labs)
    + xlab("r") 
    + ylab("G(r)") 
    + xlim(xlim)
    + ylim(ylim)
    + ggtitle(title)
  )
  print(p1)
  options(warn=0)
}



#####################################
# DATA BANKS

###
mPlot.sqa <- function(data){
  N <- length(data)
  n.x <- n.y <- 1
  if(N>=2) n.y <- 2
  if(N>=3) n.x <- 2
  par(mfrow=c(n.x, n.y), mar=c(2,4,2,1))
  for(i in 1:N)
    plot(data[[i]]$x, data[[i]]$y, t="l", ylab=paste("bank ", i, sep=" "))
	
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}


###
mPlot.results.banks <- function(fit.results, label.x="x", label.y="y", xlim=NA, ylim=NA){
  N <- length(fit.results)
  n.x <- n.y <- 1
  if(N>=2) n.y <- 2              # number of columns; maximum=2
  if(N>=3) n.x <- ceiling(N/2)   # number of rows
  
  if(is.null(dim(xlim)) || is.null(dim(ylim)))
    xlim <- ylim <- matrix(NA, nrow=N, ncol=2)  
  for(i in 1:N){
    if(any(is.na(ylim[i,])))
      ylim[i, ] <- c(min(fit.results[[i]]$curves$y-fit.results[[i]]$curves$SB), 
                     max(fit.results[[i]]$curves$y-fit.results[[i]]$curves$SB))
  }
  for(i in 1:N){
    if(any(is.na(xlim[i,])))
      xlim[i, ] <- c(min(fit.results[[i]]$x), max(fit.results[[i]]$x))
  }

  for(i in 1:N){
    fit.res <- fit.results[[i]]
    if(any(is.na(fit.res$knots$x)) || any(is.na(fit.res$knots$y))){
      fit.res$knots$x <- 0
      fit.res$knots$y <- 0  
    }
    
    if(length(fit.res$uncrt) == 1 && is.na(fit.res$pars) ){
      dat <- list(x=fit.res$x, y=fit.res$curves$y, SB=fit.res$curves$SB, 
                lambda=fit.res$fit.details$lambda, sigma=fit.res$fit.details$sigma)
      if(!is.null(fit.res$curves$corr))
        dat$y <- dat$y + fit.res$curves$corr
      fit.res$uncrt <- get.hess(dat, fit.res$knots$x, fit.res$knots$y, Gr=NA, r=seq(0, 2, 0.01), p.bkg=.5)
    }
    else{ 
      fit.res$uncrt <- list()
      fit.res$uncrt$stdev <- rep(0, length(fit.res$x))
    }
    #to avoid NOTEs in package check:
    x=y=bkg=stdev=signal=variable=value=NULL
    curves <- data.frame(x = fit.res$x, y = fit.res$curves$y-fit.res$curves$SB, 
                         bkg = fit.res$curves$bkg)				 
    melted.curves <- reshape2::melt(curves, id.vars="x")
    ribbon <- data.frame(x = fit.res$x, bkg = fit.res$curves$bkg, stdev = fit.res$uncrt$stdev,
                         signal = fit.res$curves$y-fit.res$curves$bkg)
    knots <- data.frame(x=fit.res$knots$x, y=fit.res$knots$y)
 
    assign(paste("p", i, sep=""), 
           (h <- ggplot2::ggplot(melted.curves, aes(x=x, y=value))
                + geom_line(aes(colour=variable, group=variable)) 
                + geom_ribbon(data=ribbon, aes(ymin=bkg-2*stdev, ymax=bkg+2*stdev, y=bkg), fill = "brown", alpha=0.5)
                + scale_color_manual(name="Plot legend", values = c("black", "red", "brown"))
                + geom_point(data=knots, aes(x=x, y=y), colour="red", size=2)
                + xlab(label.x) 
                + ylab(label.y) 
                + xlim(xlim[i,])
                + ylim(ylim[i,])
                + theme(legend.position = "none")))
  }

  Layout <- grid.layout(nrow = n.x+1, ncol = n.y, 
                        widths = unit(rep(2,length=n.y), rep("null", length=n.y)), 
                        heights = unit(c(rep(1,length=n.x), 0.4), rep("null", length=(n.x+1))))
     
  vplayout <- function(...) {
     grid.newpage()
     pushViewport(viewport(layout = Layout))
  }
  subplot <- function(x, y){ 
    viewport(layout.pos.row = x, layout.pos.col = y)
  }   
  vplayout()
  for(i in 1:N){
    pp <- get(paste("p", i, sep=""))
    xx <- ceiling(i/n.y)
    yy <- (i+1) %% 2 + 1
    print(pp, vp = subplot(xx, yy))
  }
  empty.curves <- reshape2::melt(data.frame(x = c(0), y = c(0), bkg = c(0))	, id.vars="x")
  options(warn=-1)
  legend <- (h <- ggplot2::ggplot(empty.curves, aes(x=x, y=value))
      + geom_line(aes(colour=variable, group=variable)) 
      + scale_color_manual(name="", values = c("black", "red"), labels = c("experimental data", "background guess (+/-2sd)"))
      + theme(legend.position=c(.5, .5), 
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks = element_blank(), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              axis.ticks.margin = unit(c(0,0,0,0), "lines"),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              line = element_blank(),
              legend.direction='horizontal',
              legend.box='vertical'              
        )
  )
  print(legend, vp = subplot(n.x+1, 1:n.y)) 
  options(warn=0)
}
