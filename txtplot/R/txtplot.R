getTicks <- function(rng){
  ## constrain pretty to pre-specified range
  ticks <- pretty(rng)
  ticks[ticks <= rng[2] & ticks >= rng[1]]
}

getRng <- function(x, perc = 0.05){
  ux <- unique(x)
  if(length(ux) == 1){
    if(ux == 0)
      rng <- c(-1,1)
    else
      rng <- ux+c(-0.4,0.4)*abs(ux)
  } else {
    rng <- range(x)
    ## extend over range (to avoid points on the axes)
    rng[1] <- rng[1]-perc*diff(rng)
    rng[2] <- rng[2]+perc*diff(rng)
  }
  rng
}

checkValid <- function(x){
  ## remove NA and NaN for plotting
  ## throw an error for infinities
  naind <- is.na(x)|is.nan(x)
  if(sum(!naind) == 0) # everything NA or NaN
    stop("nothing to plot")
  Infind <- is.infinite(x)
  if(any(Infind))
    stop("need finite values to plot")
  !naind
}

txtplot <- function(x, y = NULL, pch="*",
                    width = round(options()$width*0.8),
                    height = round(0.25*width),
                    xlab = NULL, ylab = NULL,
                    xlim = NULL, ylim = NULL){
  ## initial checks
  if(!is.null(y)){
    if(length(x) != length(y))
      stop("x and y need to have same length")
  } else {
    y <- x
    x <- 1:length(y)
  }
  if(!is.numeric(x) | !is.numeric(y))
    stop("x and y need to be of type numeric")
  indx <- checkValid(x)
  indy <- checkValid(y)
  x <- x[indy & indx]
  y <- y[indy & indx]
  ## set plotting region
  if(is.null(xlim)){
    rngx <- getRng(x)
  } else {
    if(length(xlim) != 2)
      stop("need vector of length 2 for xlim")
    if(!is.numeric(xlim))
      stop("xlim needs to be a numeric")
    ## truncate data to xlim
    ind <- x >= xlim[1] & x <= xlim[2]
    x <- x[ind]
    y <- y[ind]
    rngx <- xlim
  }
  if(is.null(ylim)){
    rngy <- getRng(y)
  } else {
    if(length(ylim) != 2)
      stop("need vector of length 2 for ylim")
    if(!is.numeric(ylim))
      stop("ylim needs to be a numeric")
    ## truncate y to ylim
    ind <- y >= ylim[1] & y <= ylim[2]
    x <- x[ind]
    y <- y[ind]
    rngy <- ylim
  }
  ## get tick marks
  xticks <- getTicks(rngx)
  yticks <- getTicks(rngy)
  ## get tick labels
  xtkch <- formatC(xticks, digits=5, width=-1)
  ytkch <- formatC(yticks, digits=5, width=-1)
  ## calculate plot margins and size of plotting area
  ## space needed for ytick marks +2 for prettyness (+ space for ylab)
  lmar <- max(nchar(ytkch)) + 2 + 2*as.numeric(!is.null(ylab))
  nch <- nchar(xtkch[length(xtkch)])
  if(is.element(nch, c(3,4))){
    rmar <- 1
  } else {
    if(nch == 5)
      rmar <- 2
    else
      rmar <- 0
  }
  pwid <- width - lmar - rmar # effective plot width
  if(!is.null(xlab)){ # x-axis label
    if(!is.character(xlab))
      stop("xlab needs to be a character vector")
    xlab <- substr(xlab, 1, pwid)
    bmar <- 3
  } else {
    bmar <- 2
  }
  phgt <- height - bmar # effective plot height
  if(!is.null(ylab)){ # y-axis label
    if(!is.character(ylab))
      stop("ylab needs to be a character vector")
    ylab <- substr(ylab, 1, phgt)
  }
  if(pwid < 6)
    stop(paste("effective plotting width <6 characters", sep=""))
  if(phgt < 6)
    stop(paste("effective plotting heigth <6 characters", sep=""))
  ## increments
  xdelta <- diff(rngx)
  ydelta <- diff(rngy)  
  ## initialize plot vector
  ch <- character((width+1)*height) # +1 for EOL marks
  ch[1:((width+1)*height)] <- " " # initialise to blank
  ## EOL marks
  ch[(1:height)*(width+1)] <- "\n" 
  ## create axes
  ind <- c(lmar:width, (phgt*(width+1)+lmar):(phgt*(width+1)+width))
  ch[ind] <- "-"
  ind <- c(lmar+(0:phgt)*(width+1), width+(0:phgt)*(width+1))
  ch[ind] <- "|"
  ind <- c(lmar, width, phgt*(width+1)+lmar, phgt*(width+1)+width)
  ch[ind] <- "+"
  ## create tick-marks
  xtck <- round((xticks-rngx[1])/xdelta*pwid)
  ind <- c(lmar+xtck, phgt*(width+1)+lmar+xtck)
  ch[ind] <- "+"
  ytck <- round((yticks-rngy[1])/ydelta*phgt)
  ind <- c(lmar+(phgt-ytck)*(width+1), width+(phgt-ytck)*(width+1))
  ch[ind] <- "+"
  ## create tick labels
  indx <- (phgt+1)*(width+1)+lmar+xtck
  xl <- strsplit(xtkch, NULL)
  for(i in 1:length(xl)){
    ln <- length(xl[[i]])
    fln <- floor(ln/2)
    z <- 1
    shift <- 0
    ## assure we do not write over plot region
    if(i == 1 & (indx[i]-fln <= (phgt+1)*(width+1))){
      shift <- indx[i]-fln-((phgt+1)*(width+1)+1)
    }
    if(i == length(xl) & (indx[i]+fln > (phgt+1)*(width+1)+width+1)){
      shift <- indx[i]+fln - ((phgt+1)*(width+1)+width)
    }
    for(j in (-fln:fln)-shift){
      if(!(ln%%2) & j == fln-shift)
        break
      ch[indx[i] + j] <- xl[[i]][z]
      z <- z+1
    }
    ## for(j in -fln:fln){
    ##   if(!(ln%%2) & j == fln)
    ##     break
    ##   ch[indx[i] + j] <- xl[[i]][z]
    ##   z <- z+1
    ## }
  }
  indy <- lmar+(phgt-ytck)*(width+1)
  yl <- strsplit(ytkch, NULL)
  for(i in 1:length(yl)){
    ln <- length(yl[[i]])
    z <- 1
    for(j in -((ln+1):2)){
      ch[indy[i] + j] <- yl[[i]][z]
      z <- z+1
    }
  }
  ## x/y-label
  if(!is.null(xlab)){
    xlabpos <- round(0.5*(rngx[2]-rngx[1])/xdelta*pwid)
    indxlab <- (phgt+2)*(width+1)+lmar+xlabpos
    xl <- strsplit(xlab, NULL)
    ln <- length(xl[[1]])
    fln <- floor(ln/2)
    z <- 1
    for(j in -fln:fln){
      if(!(ln%%2) & j == fln)
        break
      ch[indxlab + j] <- xl[[1]][z]
      z <- z+1
    }
  }
  if(!is.null(ylab)){
    ylabpos <- round(0.5*(rngy[2]-rngy[1])/ydelta*phgt)
    indylab <- 1+(phgt-ylabpos)*(width+1)
    yl <- strsplit(ylab, NULL)
    ln <- length(yl[[1]])
    fln <- floor(ln/2)
    z <- 1
    for(j in -fln:fln){
      if(!(ln%%2) & j == fln)
        break
      ch[indylab + j*(width+1)] <- yl[[1]][z]
      z <- z+1
    }
  }
  ## actual plotting
  xplt <- round((x-rngx[1])/xdelta*pwid)
  yplt <- round((y-rngy[1])/ydelta*phgt)
  ind <- (phgt-yplt)*(width+1)+lmar+xplt
  ch[ind] <- pch
  ## actual printing
  cat(ch, sep ="")
}

txtdensity <- function(x, pch = "*", width = round(options()$width*0.8),
                    height = round(0.25*width), xlab = NULL, ylab = NULL){
  dens <- density(x, from=min(x), to=max(x))
  txtplot(dens$x, dens$y, pch = pch, width = width,
          height = height, xlab = xlab, ylab = ylab)
}

txtcurve <- function(expr, from = NULL, to = NULL, n = 101,
                     pch = "*", width = round(options()$width*0.8),
                     height = round(0.25*width), xlab = NULL, ylab = NULL){
  if(any(is.null(c(from, to))))
    stop("'from' and 'to' need to be specified")
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    fcall <- paste(sexpr, "(x)")
    expr <- parse(text = fcall)
  }
  else {
    if (!((is.call(sexpr) || is.expression(sexpr)) && match("x", 
        all.vars(sexpr), nomatch = 0L))) 
      stop("'expr' must be a function, call or an expression containing 'x'")
    expr <- sexpr
  }
  x <- seq.int(from, to, length.out = n)
  y <- eval(expr, envir = list(x = x), enclos = parent.frame())
  txtplot(x, y, pch = pch, width = width, height = height, xlab = xlab, ylab = ylab)
}

txtacf <- function(x, pch = "*", lag.max = 20, type = c("correlation", "covariance", "partial"),
                   na.action = na.fail, demean = TRUE, 
                   width = round(options()$width*0.8), height = round(0.25*width),
                   xlab = NULL, ylab = NULL){
  acfplot <- acf(x, lag.max = lag.max, type = type, na.action = na.action,
                 demean = demean, plot = FALSE)
  func <- function(x,height){
    if(x > 0)
      out <- seq(0, x, length = height)
    if(x < 0)
      out <- seq(x, 0, length = height)
    out
  }
  lst <- lapply(acfplot$acf, func, height=height)
  yy <- do.call("c", lst)
  xx <- rep(acfplot$lag, each = height)
  txtplot(xx, yy, width = width, height = height, pch = pch, xlab = xlab, ylab = ylab)
}

txtbarchart <- function(x, pch = "*",
                        width = round(options()$width*0.8), height = round(0.25*width),
                        ylab = NULL){
  if(!is.factor(x))
    stop("x needs to be a factor")
  sm <- summary(x)
  xx0 <- 1:length(sm)
  yy0 <- sm/sum(sm)*100
  func <- function(x,height){
    seq(0, x, length = height)
  }
  lst <- lapply(yy0, func, height=height)
  yy <- do.call("c", lst)
  xx <- rep(xx0, each = height)
  txtplot(xx, yy, pch=pch, width = width, height = height-1, ylab = ylab)
  nams <- attr(sm, "names")
  drawLegend(nams, width)
}

txtboxplot <- function(..., range = 1.5, legend = NULL, xlab = NULL,
                       width = round(options()$width*0.8)){
  ## preliminary steps
  objs <- list(...)
  nobj <- length(objs)
  nams <- as.character(match.call())[2:(nobj + 1)]
  ## write everything in one vector
  boxlist <- lapply(objs, function(x) boxplot.stats(x, coef = range)$stats)
  x <- do.call("c", boxlist)
  ## determine whether to use legend
  if(is.null(legend))
    legend <- ifelse(nobj == 1, FALSE, TRUE)
  ## determine plotting region
  rng <- getRng(x, perc=0.08)
  delta <- diff(rng)
  ## get tick marks
  xticks <- getTicks(rng)
  ## get tick labels
  xtkch <- formatC(xticks, digits=5, width=-1)
  ## initialize tick labels
  ch0 <- ch1 <- ch2 <- character(width+1)
  ch0[1:(width+1)] <- ch1[1:(width+1)] <- ch2[1:(width+1)] <- " "
  ## EOL marks
  ch0[width+1] <- ch1[width+1] <- ch2[width+1] <- "\n"
  ## xlabels
  if(!is.null(xlab)){
    xlabpos <- round(0.5*(rng[2]-rng[1])/delta*width)
    indxlab <- xlabpos
    ## truncate to plotting region
    xlab <- substr(xlab, 1, width) 
    xl <- strsplit(xlab, NULL)
    ln <- length(xl[[1]])
    fln <- floor(ln/2)
    z <- 1
    for(j in -fln:fln){
      if(!(ln%%2) & j == fln)
        break
      ch0[xlabpos + j] <- xl[[1]][z]
      z <- z+1
    }
    cat(ch0, sep="")
  }
  ## create axis
  ind <- c(2:width)
  ch2[ind] <- "-"
  ## create tick-marks
  xtck <- round((xticks-rng[1])/delta*(width-1))
  ind <- c(1+xtck)
  ch2[ind] <- "+"
  ind <- c(2, width)
  ch2[ind] <- "|"
  ## create tick labels
  indx <- 1+xtck
  xl <- strsplit(xtkch, NULL)
  for(i in 1:length(xl)){
    ln <- length(xl[[i]])
    fln <- floor(ln/2)
    z <- 1
    shift <- 0
    ## assure we do not write over plot region
    if(i == 1 & (indx[i]-fln < 0)){
      shift <- indx[i]-fln-1
    }
    if(i == length(xl) & (indx[i]+fln > width)){
      shift <- indx[i]+fln - width
    }
    for(j in (-fln:fln)-shift){
      if(!(ln%%2) & j == fln-shift)
        break
      ch1[indx[i] + j] <- xl[[i]][z]
      z <- z+1
    }
  }
  cat(ch1, sep="")
  cat(ch2, sep="")  
  
  ## produce boxplots
  for(k in 1:nobj){
    if(!is.numeric(objs[[k]]))
      stop("vectors to be plotted need to be of type numeric")
    bstats <- boxplot.stats(objs[[k]], coef = range)$stats
    if(legend)
      k0 <- k
    else 
      k0 <- NULL
    boxcore(bstats, width, rng, delta, k0)
  }
  if(!is.null(k0)){
    drawLegend(nams, width)
  }
}

boxcore <- function(bstats, width, rng, delta, no = NULL){
  inc <- 1
  ch <- character(3*(width+1))
  ch[1:(3*(width+1))] <- " "
  ch[(1:3)*(width+1)] <- "\n"
  marks <- round((bstats-rng[1])/delta*(width-1))
  ## first line
  ch[inc+marks[2]:marks[4]] <- "-"
  ch[inc+marks[2:4]] <- "+"
  ## middle line
  if(!is.null(no))
    ch[1 + width + 1] <- no
  ch[inc+marks[1]:marks[2]+width+1] <- "-"
  ch[inc+marks[4]:marks[5]+width+1] <- "-"
  ch[inc+marks[2:4]+width+1] <- "|"
  ## last line
  ch[inc+marks[2]:marks[4]+(width+1)*2] <- "-"
  ch[inc+marks[2:4]+(width+1)*2] <- "+"
  cat(ch, sep="")
}

insEOL <- function(x, width){
  ## insert line endings into a character vector
  lenx <- length(x)
  lenins <- floor(lenx/width)
  if(lenins > 0)
    out <- character(lenx + lenins+1)
  else
    return(c(x, "\n"))
  ind <- setdiff(1:(length(out)-1), (1:lenins)*(width+1))
  out[ind] <- x
  ind <- c((1:lenins)*(width+1), lenx + lenins+1)
  out[ind] <- "\n"
  out
}

drawLegend <- function(nams, width){

  leg <- paste(1:length(nams), nams, sep="=")
  leg <- paste(leg, collapse = ", ")
  leg <- insEOL(strsplit(leg, NULL)[[1]], width)
  cat("Legend: ")
  if(length(leg)>width)
    cat("\n")
  cat(leg, sep="")
}

