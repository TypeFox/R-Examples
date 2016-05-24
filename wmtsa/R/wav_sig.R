################################################
## WMTSA series functionality
##
##  Functions:
##
##    create.signalSeries
##    make.signal
##
################################################

###
# create.signalSeries
##

"create.signalSeries" <- function(x=NULL, position=list(from=1,by=1,units=character()),
  from=NULL, by=NULL, to=NULL, length.out=NULL,
  units=character(), title.data=character(), documentation=character(), na.rm=TRUE)
{
  if (is(x,"signalSeries"))
    return(x)
  if (is(x,"ts")){
    title.data <- deparseText(substitute(x))
    position   <- list(from=start(x)[1], by=deltat(x), units=character())
    x          <- as.vector(x)
    length.out <- length(x)
  }

  # define local functions
  "fill.missing.position" <- function(x)
  {
    defaults <- list(from=1, by=1, units=character(0))
    matches  <- which(as.logical(charmatch(names(defaults), names(x))))

    if (length(matches) < 1) stop("Invalid position name")
    else position <- c(defaults[seq(along=defaults)[-matches]], x)

    position
  }

  # check for blank data set. if blank, return the
  # object without any position or amplitude data,
  # but put in place the units, title, and documentation
  if (is.null(x)){
    y                <- signalSeries()
    y@units.position <- position$units
    y@units          <- units
    y@title          <- title.data
    y@documentation  <- documentation

    return(y)
  }

  # check for NA

  na.indices <- is.na(x)

  if (any(na.indices)){
    if(!na.rm)
      stop("NaNs not allowed in x. Try setting: na.rm=TRUE")
    x <- x[!na.indices]
    warning("NaNs removed from x for conversion to signalSeries class")
  }

  # check for NULL arguments
  from.is.null   <- is.null(from)
  to.is.null     <- is.null(to)
  length.is.null <- is.null(length.out)
  by.is.null     <- is.null(by)

  # check whether from is an index vector of
  # the form a:b
  is.from.an.index.vector <- (is.integer(from) & isVectorAtomic(from) & length(from) > 1)

  # verify appropriate ranges
  if (!length.is.null){
    if (length.out < 1)
      stop("Length argument must be a positive integer")
  }

  # fill in missing position objects
  position <- fill.missing.position(position)

  # set from if not specified
  if (from.is.null)
    from <- position$from

  # obtain x class and convert
  # to a one dimensional vector
  # (not necessarily of numeric class,
  # it can be ts class)
  cls <- class(x)

  if (cls == "signalSeries"){
    x@title <- title.data
    x@documentation <- documentation
    return(x)
  }

  if (isVectorAtomic(x))
    y <- as.vector(x)
  else
    stop("Input must be a vector, a single-column or single-row matrix,
      or a named vector")

  if(is(x, "ts")){

    if (is.from.an.index.vector){
      # obtain time series specifications
      tsp.data   <- tsp(y)
      start.data <- tsp.data[1]
      end.data   <- tsp.data[2]
      freq.data  <- tsp.data[3]

      # perform linear interpolation to calculate
      # dates corresponding to the index vector
      date.range <- approx(x=c(1,length(y)), y=c(start.data, end.data),
        xout=c(from[1],from[length(from)]), method="linear")

      # create a time series beginning at same time as ship
      y <- ts(data=y[from], start=date.range$y[1],
               frequency=freq.data)
    }

    #y <- ts.update(y)

  }
  else{

    if (is.from.an.index.vector)
      index <- from[from <= length(y)]
    else{
      # calculate length of vectorized data
      ly <- length(y)

      # develop position vector
      pos.all <- seq(from=position$from, by=position$by, length=ly)
      if (!to.is.null && by.is.null){
	    if (from < to)
	      index <- which(pos.all >= from & pos.all <= to)
	    else
	      index <- which(pos.all >= to & pos.all <= from)

        if(length(index) > length.out)
	      index <- index[seq(length.out)]
      }
      else{
	    if (to.is.null)
	      to <- pos.all[ly]
	    if (by.is.null)
	      by <- sign(to-from)*position$by

        pos.larger   <- pos.all[pos.all >= from]
	    closest.from <- pos.larger[order(abs(pos.larger - from))[1]]

        if (by.is.null){
	      if (!anyMissing(closest.from)){
	        if (closest.from != from)
	          extract <- seq(from=closest.from, by=by, to=to)
	        else
	          extract <- seq(from=from, by=by, to=to)
	      }
	      else
	        extract <- NA
	    }
	    else
	      extract <- seq(from=from, by=by, to=to)

	    if (from <= max(pos.all)){
	      index <- pmatch(extract, pos.all)
	      index <- index[!is.na(index)]
	      index <- index[seq(length=min(length.out,length(index)))]
	    }
	    else
	      index <- NA
      }
    }

    if (anyMissing(index) || any(is.null(index)))
      stop("Indexing out of range")

    y <- y[index]

    # finally, convert to a signalSeries class
    new.position=(index-1)*position$by + position$from
    y <- signalSeries(data=y,
		      from=new.position[1], by=diff(new.position[1:2]),
		      units=units, units.position=position$units)
  }

  # add title and documentation slots
  if (length(y@units) < 1 && !is.null(units))
    y@units <- units
  if (length(y@title) < 1 && !is.null(title.data))
    y@title <- title.data
  if (length(y@documentation) < 1 && !is.null(documentation))
    y@documentation <- documentation
  if (length(y@units.position) < 1 && !is.null(position$units))
    y@units.position <- position$units

  y
}

###
# make.signal
##

"make.signal" <- function(name, n=1024, snr=Inf)
{

 ".wave.demo.signals" <- c("dirac", "kronecker", "heavisine", "bumps", "blocks",
  "doppler", "ramp", "cusp", "crease", "sing", "hisine",
  "losine", "linchirp", "twochirp", "quadchirp",
  "mishmash1", "mishmash2", "mishmash3", "levelshift",
  "jumpsine", "gauss", "patches",
  "linear", "quadratic", "cubic")

  x <- (0:(n-1.))/n
  z <- switch(name,
   dirac=n*(x == floor(.37*n)/n),
   kronecker=(x == floor(.37*n)/n),
   heavisine=4*sin(4*pi*x)-sign(x-.3)-sign(.72-x),
   bumps={
     pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
     hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
     wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
     y <- rep(0, n)
     for(j in 1:length(pos)) y <- y+hgt[j]/(1+abs((x-pos[j]))/wth[j])^4
     y
   },
   blocks={
     pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
     hgt <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1,2.1, -4.2)
     y <- rep(0, n)
     for(j in 1:length(pos)) y <- y+(1+sign(x-pos[j]))*hgt[j]/2
     y
   },
   doppler=sqrt(x*(1-x))*sin((2*pi*1.05)/(x+.05)),
   ramp=x-(x >= .37),
   cusp=sqrt(abs(x-.37)),
   crease=exp(-4*abs(x-.5)),
   sing=1/abs(x-(floor(n*.37)+.5)/n),
   hisine=sin(pi*n*.6902*x),
   midsine=sin(pi*n*.3333*x),
   losine=sin(pi*n*.03*x),
   linchirp=sin(.125*pi*n*x^2),
   twochirp=sin(pi*n*x^2) + sin((pi/3)*n*x^2),
   quadchirp=sin((pi/3)*n*x^3),
   # QuadChirp + LinChirp + HiSine
   mishmash1=sin((pi/3)*n*x^3) + sin(pi*n*.6902*x) + sin(pi*n*.125*x^2),
   # QuadChirp + LinChirp + HiSine + Bumps
   mishmash2={		# wernersorrows
     y   <- sin(pi*(n/2)*x^3)+sin(pi*n*.6902*x)+sin(pi*n*x^2)
     pos <- c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
     hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
     wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
     for(j in 1:length(pos)) y <- y + hgt[j]/(1+abs((x-pos[j])/wth[j]))^4
     y
   },
   # QuadChirp + MidSine + LoSine + Sing/200.
   mishmash3=sin((pi/3)*n*x^3) + sin(pi*n*.3333*x) + sin(pi*n*.03*x) +
        (1/abs(x-(floor(n*.37)+.5)/n))/(200.*n/512.),
   gauss=dnorm(x, .3, .025),
   jumpsine=10.*(sin(4*pi*x) + as.numeric(x >= 0.625 & x < 0.875)),
   levelshift=as.numeric(x >= 0.25 & x < 0.39),
   patches={
     if(n<16) stop("n must be >= 16 to generate patches\n")
     J <- ilogb(n, base=2)
     y <- rep(0., n)
     for(j in 0:(J-4.)) y[(1:2^j)+3.*2.^(j+2.)] <- 1.
     y
   },
   linear=2.*x-1.,
   quadratic=4.*(1.-x)*x,
   cubic=64.*x*(x-1.)*(x-.5)/3.,
   stop("Unknown signal name.  Allowable names are:\n",
	paste(.wave.demo.signals, collapse=", ")))

   if (snr > 0)
     z <- z + rnorm(n)*sqrt(var(z))/snr

   z <- signalSeries(data=z, from=0.0, by=1.0/n)
   z@title <- name
   z
}
