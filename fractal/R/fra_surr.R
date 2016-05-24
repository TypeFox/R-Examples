################################################
## FRACTAL surrogate data constructor
## functions and corresponding methods
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class surrogate
## Constructor function: surrogate
## Methods:
##
##   eda.plot.surrogate
##   plot.surrogate
##   print.surrogate
##
################################################

###
# surrogate
###

"surrogate" <- function(x, method="ce", sdf=NULL, seed=0)
{
  # obtain series name
  if (is(x,"signalSeries"))
    series.name <- ifelse1(length(x@title), x@title, deparseText(substitute(x)))
  else
    series.name <- deparseText(substitute(x))

  series   <- asVector(x)
  n.sample <- length(x)

  # check seed
  checkScalarType(seed,"integer")
  seed <- as.integer(abs(seed))

  # map the surrogate method
  checkScalarType(method,"character")
  method <- match.arg(lowerCase(method), c("phase","aaft","dh","ce"))

  if (method == "ce"){

    # calculate the spectral density function
    if (is.null(sdf)){

      h   <- sapa::taper(type="sine", n.taper=min(c(5,n.sample)), n.sample=n.sample, normalize=TRUE)
      sdf <- Re(sapa::SDF(x, recenter=TRUE, taper.=h, method="multitaper"))
    }

    if (!is(sdf,"SDF"))
      stop("SDF must be of class \"SDF\"")
    if (!attr(sdf, "single.sided"))
      stop("SDF must be single sided")
    if (1 / (2 * attr(sdf,"deltaf") * attr(sdf,"deltat")) != attr(sdf,"n.sample"))
      stop("\n\n\tFor the circulant embedding technique, the input SDF must\n",
       "\thave a frequency resolution of 1/(2N) where N is the\n",
       "\tnumber of samples in the original time series. The input SDF\n",
       "\tcurrently contains a spectral resolution of ", attr(sdf,"deltaf"),
       "\n\tand it should be ", 1 / (2*attr(sdf,"n.sample")* attr(sdf,"deltat")), ". Please recreate\n",
       "\tthe SDF ala SDF(x, npad=", 2*attr(sdf,"n.sample"), ", ...) to accomplish this\n",
       "\tand try again.\n")

    # make sure to multiply the SDF estimate by
    # the sampling interval since the C-code
    # assumes a unit sampling interval
    S <- matrix(as.vector(sdf)) / attr(sdf, "deltat")
    z <- as.vector(itCall("RS_fractal_bootstrap_circulant_embedding",
      S, as.integer(seed)))
      #COPY=rep(FALSE,2),
      #CLASSES=c("matrix","integer"),
      #PACKAGE="ifultools"))

    # restore mean if recentered
    if (attr(sdf, "recenter"))
      z <- z + mean(series)
  }
  else if (method == "dh"){

    z <- as.vector(itCall("RS_fractal_bootstrap_davison_hinkley",
      as.numeric(series), as.integer(seed)))
      #COPY=rep(FALSE,2), #CLASSES=c("matrix","integer"),
      #PACKAGE="ifultools"))
  }
  else{

    z <- as.vector(itCall("RS_fractal_bootstrap_theiler",
      as.numeric(series), as.integer(ifelse1(method == "phase", 0, 1)), as.integer(seed)))
      #COPY=rep(FALSE,3), #CLASSES=c("matrix", "integer","integer"),
      #PACKAGE="ifultools"))
  }

  # assign class
  oldClass(z) <- "surrogate"

  attr(z, "method")      <- method
  attr(z, "sub.method")  <- ifelse1(method=="aaft", "phase", NULL)
  attr(z, "n.sample")    <- n.sample
  attr(z, "series.name") <- series.name
  attr(z, "series")      <- series
  attr(z, "seed")        <- ifelse1(seed > 0, as.integer(seed), NULL)

  z
}

###
# eda.plot.surrogate
###

"eda.plot.surrogate" <- function(x,
  show.="surrogate", col.series="black", col.surrogate="red", cex=NULL, ...)
{
  xatt <- attributes(x)
  series.name <- xatt$series.name
  show. <- match.arg(show., c("series","surrogate","both"))

  if (show. == "both"){

    nr <- 4
    nc <- 2
    if (is.null(cex))
      cex <- 0.4
  }
  else{
  	nr <- 2
  	nc <- 2
    if (is.null(cex))
      cex <- 0.85
  }

  old.plt <- splitplot(nr,nc,1)
  on.exit(par(old.plt))

  types <- c("time","pdf","sdf","lag")

  for (i in seq(along=types)){

    iplot <- ifelse1(show.=="both", 2*i-1, i)

    if (i > 1)
      splitplot(nr,nc,iplot)

    plot(x, show.=show., add=TRUE,
      type=types[i], cex=cex, stack=FALSE,
      split.=c(nr, nc, iplot),
      col.series=col.series, col.surrogate=col.surrogate, ...)
  }

  invisible(NULL)
}

###
# plot.surrogate
###

"plot.surrogate" <- function(x,
  show.="surrogate", type="time", stack=TRUE,
  xlab=NULL, ylab=NULL,
  add=FALSE, cex=1,
  adj.main=1, line.main=0.5, split.=NULL,
  col.series="black", col.surrogate="red", ...)
{
  if (!is.null(split.) && (!is.vector(split.) || length(split.) != 3))
    stop("split must be a 3 element vector")

  type  <- match.arg(type, c("time","pdf","sdf","lag"))
  show. <- match.arg(show., c("series","surrogate","both"))

  # cannot use stackPlot for density and lag plots since the
  # x-axis differs in each case for the original and surrogate series
  if (show. == "both" && is.element(type, c("pdf","lag")))
    stack <- FALSE

 	# split the plotting region if not done so exogenously
 	if (is.null(split.)){
 	  nr      <- ifelse1(show. == "both", ifelse1(stack, 1, 2), 1)
	  old.plt <- splitplot(nr, 1, 1)
          on.exit(par(old.plt))
 	}

  # initialize variables
  xatt        <- attributes(x)
  series.name <- xatt$series.name
  surrname    <- paste(upperCase(xatt$method),"Surrogate")
  series      <- xatt$series
  surrog      <- asVector(x)
  series.time <- ifelse1(is(series, "signalSeries"),
    as(positions(series),"numeric"), time(series))

  # plot data
  if (type == "time"){

    if (is.null(xlab))
      xlab <- "Time"
    if (is.null(ylab))
      ylab <- series.name

    if (show. == "both"){

      if (stack){

        stackPlot(x=as.vector(series.time), data.frame(series, surrog),
          xlab=xlab, ylab=list(text=c(series.name, surrname),
          col=c(col.series, col.surrogate)), add=add, ...)
        title(main="Time History", adj=adj.main, line=line.main)
      }
      else{

        plot(x=series.time, y=series, col=col.series, type="l", cex=cex,
          xlab=xlab, ylab=ylab, ...)
        title(main="Time History", adj=adj.main, line=line.main)

        if (!is.null(split.))
          splitplot(split.[1],split.[2],split.[3]+1)
        else
          splitplot(nr, 1, 2)

       	plot(x=series.time, y=surrog, col=col.surrogate, type="l", cex=cex,
           xlab=xlab, ylab=surrname, ...)
      }
    }
    else{
     	plot(x=series.time, y=ifelse1(show.=="surrogate",surrog, series),
     	  col=ifelse1(show.=="surrogate",col.surrogate, col.series), type="l", cex=cex,
        xlab=xlab, ylab=ifelse1(show.=="surrogate", surrname, series.name), ...)
      title(main="Time History", adj=adj.main, line=line.main)
    }

  }
  else if (type == "pdf"){

    p.series <- density(series, ...)
    p.surrog <- density(surrog, ...)
    main <- "PDF"

    if (show. == "both"){

     	plot(p.series, main="", ylab=paste(series.name, "Density"), col=col.series)
      title(main=main, adj=adj.main, line=line.main)

      if (!is.null(split.))
        splitplot(split.[1],split.[2],split.[3]+1)
      else
        splitplot(nr, 1, 2)

      plot(p.surrog, main="", ylab=paste(surrname, "Density"), col=col.surrogate, ...)
    }
    else{

     	plot(ifelse1(show.=="surrogate", p.surrog, p.series), main="",
     	  ylab=paste(ifelse1(show.=="surrogate", surrname, series.name), "Density"),
     	  col=ifelse1(show.=="surrogate",col.surrogate, col.series), ...)
      title(main=main, adj=adj.main, line=line.main)
    }

  }
  else if (type == "sdf"){

    if (is.null(xlab))
      xlab <- "Frequency (Hz)"

    main <- "SDF"

    S.series  <- sapa::SDF(series, method="multitaper")
    S.surrog  <- sapa::SDF(surrog, method="multitaper")
    frequency <- as.vector(attr(S.series,"frequency"))

    if (show. == "both"){

      if (stack){

        stackPlot(x=frequency, y=decibel(data.frame(as.numeric(S.series), as.numeric(S.surrog))),
          col=c(col.series, col.surrogate), xlab=xlab,
          ylab=list(text=paste(c(series.name, surrname), "SDF (dB)", sep="\n"),
            col=c(col.series, col.surrogate)),
          add=add, cex=cex, ...)
        title(main=main, adj=adj.main, line=line.main)
      }
      else{

        plot(S.series, col=col.series, cex=cex, add=add, ylab=paste(series.name, "MT SDF"), ...)
        title(main=main, adj=adj.main, line=line.main)

        if (!is.null(split.))
          splitplot(split.[1],split.[2],split.[3]+1)
        else
          splitplot(nr, 1, 2)

        plot(S.surrog, col=col.surrogate, cex=cex, add=add, ylab=paste(surrname, "MT SDF"), ...)
      }
    }
    else{

     	plot(ifelse1(show.=="surrogate",S.surrog, S.series),
     	  ylab=paste(ifelse1(show.=="surrogate", surrname, series.name), "MT SDF"),
     	  col=ifelse1(show.=="surrogate", col.surrogate, col.series), ...)
      title(main=main, adj=adj.main, line=line.main)
    }

  }
  else{ # Lag plot

    emb.series <- embedSeries(series, ...)
    emb.surrog <- embedSeries(surrog, ...)

    if (show. == "both"){

     	plot(emb.series, col=col.series, ...)
      title(main="Lag Plot", adj=adj.main, line=line.main)

      if (!is.null(split.))
        splitplot(split.[1],split.[2],split.[3]+1)
      else
        splitplot(nr, 1, 2)

      plot(emb.surrog, col=col.surrogate, ...)
    }
    else{

     	plot(ifelse1(show.=="surrogate", emb.surrog, emb.series),
     	  col=ifelse1(show.=="surrogate", col.surrogate, col.series), ...)
      title(main="Lag Plot", adj=adj.main, line=line.main)
    }
  }

  invisible(NULL)
}

###
# print.surrogate
###

"print.surrogate" <- function(x, justify="left", sep=":", ...)
{
  xatt   <- attributes(x)
  method <- lowerCase(xatt$method)

  method <- switch(method,
    phase = "Random Phase",
    aaft  = "Amplitude Adjusted Fourier Transform",
    ce    = "Circulant Embeddding",
    dh    = "Davison-Hinkely Random Phase and Amplitude",
      "other")

  main <- paste("Surrogate data generation for", xatt$series.name)

  z <- list(
    "Method"=method,
    "Sub-method"=ifelse1(xatt$method == "aaft", "Random Phase", NULL),
    "Length of series"=xatt$n.sample,
    "Seed"=xatt$seed)

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}








