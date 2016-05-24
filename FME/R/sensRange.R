## -----------------------------------------------------------------------------
## Sensitivity ranges
## -----------------------------------------------------------------------------

sensRange <- function( func, parms = NULL, sensvar = NULL, dist = "unif",
         parInput = NULL, parRange = NULL, parMean = NULL, parCovar = NULL,
         map = 1, num = 100, ...) {

  vec2mat <- function(vec) { # make a matrix of a vector
    NN <- names(vec)
    mat <- matrix(data = vec, nrow = 1)
    colnames(mat) <- NN
    mat
  }
  if(is.vector(parInput))
    parInput <- vec2mat(parInput)
  if(is.vector(parRange))
    parRange <- vec2mat(parRange)
  if(is.vector(parMean))
    parMean  <- vec2mat(parMean)
  if(is.vector(parCovar))
    parCovar <- vec2mat(parCovar)

  if(is.vector(parInput))
    parInput <- matrix(data = parInput, nrow = 1)
  if (is.null(parms) & ! is.null(parInput))
    parms <- parInput[1,]
  if (is.null(parms))
    parms <- parMean
  if (is.null(parms)){
    if (is.vector(parRange))
      parms <- mean(parRange)
    else parms <- rowMeans(parRange)
  }
  if (is.null(parms))
    stop ("'parms' not known")
  if (is.matrix(parms) && nrow(parms) > 1)
    stop ("'parms' should be a vector")

  Solve <- function(parms) func(parms, ...)

  ip <- NULL
  if (! is.null(parInput)) {
    dist <- "input"
    nr <- nrow(parInput)
    num <- min(num, nr)
    if (num == nr)
       ip <- 1: nr
    else ip <- sample(1:nr, size = num, replace = FALSE)
    parset <-as.matrix(parInput[ip,])
    if(is.null(parms))
      parms <- parInput[1,]
  }

  Parms <- parms

  ## reference run
  yRef  <- Solve(Parms)

  if (is.vector(yRef)) {
    ynames <- names(yRef)
    yRef <- matrix(data = yRef, nrow = 1)
    colnames(yRef) <- ynames
  } else
  if (is.data.frame(yRef))
    yRef <- as.matrix(yRef)

  ## check sensitivity variables
  if (is.null(sensvar)) {
    ivar    <- 1:ncol(yRef)
    if (! is.null(map))
      ivar <- ivar[-map]
    sensvar <- colnames(yRef)[ivar]
    if(is.null(sensvar))
      sensvar <- ivar
  } else {
    ivar  <- findvar(yRef[1,], sensvar, "variables")
    if (! is.character(sensvar)) {# try to create names rather than nrs
      sv <- sensvar
      sensvar<-colnames(yRef)[ivar]
      if (is.null(sensvar))
        sensvar<-sv
    }
  }
  if (is.null(map))
    map   <- 1:nrow(yRef)
  else map <- yRef[,map]

  nout  <- length(ivar)
  if (nout == 0)
    stop (" should select at least ONE output variable - set map = NULL")
  ndim  <- nrow(yRef)
  grvar <- expand.grid(map, sensvar)
  if (ndim ==1)
    svar <- sensvar
  else svar <- paste(grvar[, 2], grvar[, 1], sep="")

  YREF  <- as.vector(yRef[, ivar])
  Sens  <- matrix(data=NA, nrow = num, ncol = length(YREF))

  ## sensitivity parameters
  senspar <- NULL
  senspar <- colnames(parMean)
  if (is.null(senspar)) senspar <- rownames(parRange)
  if (is.null(senspar)) senspar <- rownames(parCovar)
  if (is.null(senspar)) senspar <- colnames(parCovar)
  if (is.null(senspar)) senspar <- colnames(parInput)
  if (is.null(senspar)) senspar <- names(Parms)
  if (is.null(senspar)) senspar <- 1:length(Parms)

  npar  <- length(senspar)

  ipar <- findvar(Parms, senspar, "parameters")
  pp   <- unlist(Parms)[ipar]

  ## sanity checks for random parameters
  if (dist == "norm" && (is.null(parMean) | is.null(parCovar)))
    stop("parMean and parCovar should be given if dist = norm")
  if(!(dist  %in% c("norm", "input")) && is.null(parRange))
    stop("parRange should be given if dist = unif, grid or latin")

  ## generate random parameters
  if (dist == "norm")
    parset <- Norm (parCovar=parCovar, parMean=parMean, parRange, num) else
  if (dist == "unif")
    parset <- Unif(parRange, num)                      else
  if (dist == "latin")
    parset <- Latinhyper(parRange, num)                else
  if (dist == "grid")
    parset <- Grid(parRange, num)                      else
  if (dist != "input" )
    stop("dist should be one of 'norm', 'unif', 'latin', 'grid' or 'input'")


  ## The sensitivity output
  colnames(Sens) <- svar

  for (i in 1:num) {
    if (prod(Parms[ipar] == parset[i, ]) == 0) { # no need to run model again if same parameter value (e.g. MCMCrun!)
      Parms[ipar]  <- parset[i,]
      yRef <- Solve(Parms)
    }
    if (is.vector(yRef))
      Sens[i, ] <- yRef[ivar]
    else Sens[i,] <- as.vector(unlist(yRef[,ivar]))   # unlist in case it is a data.frame
  }
  sens<- data.frame(cbind(parset, Sens))
  class(sens) <- c("sensRange", "data.frame")
  attr(sens, "pset") <- ip   # if parInput: which parameters were drawn...
  attr(sens, "npar") <- ncol(parset)
  attr(sens, "x")    <- map
  attr(sens, "nx")   <- length(map)
  attr(sens, "var")  <- sensvar
  return (sens)
}

## -----------------------------------------------------------------------------
## S3 methods of sensRange
## -----------------------------------------------------------------------------

summary.sensRange<-function(object, ...) {

  npar <- attr(object, "npar")
  sens <- as.matrix(object[, -(1:npar)])
  x    <- attr(object, "x")
  names(x) <- NULL
  nx   <- attr(object, "nx")
  varnames  <- attr(object, "var")

  if (ncol(sens)>1)
    SumSens <- data.frame(
      x    = x,
      Mean = apply(sens, 2, FUN = mean),
      Sd   = apply(sens, 2, FUN = sd),
      Min  = apply(sens, 2, FUN = min),
      Max  = apply(sens, 2, FUN = max),
      q05  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.05)),
      q25  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.25)),
      q50  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.5)),
      q75  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.75)),
      q95  = apply(sens, 2, FUN=function(x)quantile(x, probs = 0.95))
  ) else
    SumSens <- data.frame(
      x    = x,
      Mean = mean(sens, na.rm = TRUE),
      Sd   = NA,
      Min  = min(sens, na.rm = TRUE),
      Max  = max(sens, na.rm = TRUE),
      q05  = quantile(sens, probs = 0.05, na.rm = TRUE),
      q25  = quantile(sens, probs = 0.25, na.rm = TRUE),
      q50  = quantile(sens, probs = 0.5, na.rm = TRUE),
      q75  = quantile(sens, probs = 0.75, na.rm = TRUE),
      q95  = quantile(sens, probs = 0.95, na.rm = TRUE)
    )

  rownames(SumSens) <- colnames(sens)
  attr(SumSens, "var") <- attr(object, "var")
  attr(SumSens, "nx")  <- attr(object, "nx")
  class(SumSens)<-c("summary.sensRange", "data.frame")

  return(SumSens)
}

## -----------------------------------------------------------------------------

plot.sensRange<-function(x, xyswap = FALSE, which = NULL,
  ask = NULL, ...) {
  npar <- attr(x, "npar")
  nx   <- attr(x, "nx")
  varnames  <- attr(x, "var")
  X    <- attr(x, "x")
  sens <- x[, -(1:npar)]

  dots   <- list(...)
  nmdots <- names(dots)
  Main   <- is.null(dots$main)
  dots$type <- if(is.null(dots$type)) "l" else dots$type
  Ylim   <- is.null(dots$ylim)
  Xlim   <- is.null(dots$xlim)

  Select <- selectvar(which, varnames, Nall = TRUE)

  if (nx > 1) {
    ## Set par mfrow and ask.
    ask <- setplotpar (nmdots, dots, length(Select), ask)

    ## interactively wait if there are remaining figures
    if (ask) {
        oask <- devAskNewPage(TRUE)
 	      on.exit(devAskNewPage(oask))
    }
    for (i in Select){
      ii <- ((i-1)*nx+1):(i*nx)
      if (!xyswap) {
        do.call("matplot", c(alist(X, t(sens[, ii])), dots))
      } else {
        if (Ylim) dots$ylim <- rev(range(X))
        do.call("matplot", c(alist(t(sens[, ii]), X), dots))
      }
    }  # end for
  } else
    boxplot(sens[, Select], ...)
}

## -----------------------------------------------------------------------------

plot.summary.sensRange<-function(x, xyswap = FALSE, which = NULL,
   legpos = "topleft", col = c(grey(0.8), grey(0.7)), quant = FALSE,
   ask = NULL, obs = NULL, obspar = list(), ...) {

  nx    <- attr(x, "nx")
  varnames   <- attr(x, "var")

## check observed data - can be a list
    nobs <- 0

    if (! is.null(obs)) {

      if (!is.data.frame(obs) & is.list(obs)) { # a list with different data sets
       Obs <- obs
       obs <- Obs[[1]]  
       obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
       if (length(Obs) > 1)
         for ( i in 2 : length(Obs)) {
           obs <- mergeObs(obs, Obs[[i]])
           obs.pos <- rbind(obs.pos, c(obs.pos[nrow(obs.pos),2] +1, nrow(obs)))
         }
       obsname <- colnames(obs) 
      } else {
       if (is.character(obs[,1]) | is.factor(obs[,1]))   # long format - convert
          obs <- convert2wide(obs)
       obsname <- colnames(obs) 
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
       obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
      }                       
    DD <- duplicated(obsname)
    if (sum(DD) > 0)  
      obs <- mergeObs(obs[,!DD], cbind(obs[,1],obs[,DD]))
    nobs <- nrow(obs.pos)   
    }

    dots   <- list(...)
    nmdots <- names(dots)
    Which <- which
    if (is.null (Which) & ! is.null(obs))
      Which <- which (varnames %in% obsname)
    
    Select <- selectvar(Which, varnames, Nall = TRUE)
    if (length(Select) == 0)
      stop("There is nothing to plot")
      
    np <- length (Select)
  
## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")    
    
    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot

    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab, "x"      , np)
    dotmain$ylab <- expanddots(dots$ylab, "y"      , np)
    dotmain$main <- expanddots(dots$main, varnames[Select], np)

    # ylim and xlim can be lists and are at least two values
    yylim  <- expanddotslist(dots$ylim, np)
    xxlim  <- expanddotslist(dots$xlim, np)

    # point/line parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]

## Position of variables in "obs" (NA = not observed)
    if (nobs > 0) {
      ObsWhich <- selectvar(varnames[Select], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA  
      Obspar <- setdots(obspar, nobs)
    } else 
      ObsWhich <- rep(NA, np)

  if (nx > 1)  {    # summary of a times series or a profile...
    ## Set par mfrow and ask.
    ask <- setplotpar (nmdots, dots, length(Select), ask)

    ## interactively wait if there are remaining figures
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }


    for (ij in 1:np){
      i <- Select[ij]
      io <- ObsWhich[ij]   # position of variable in 'obs'

      ii <- ((i-1)*nx+1):(i*nx)
      X <- x[ii,]
      if (quant) {
        xmin <- X$q05
        xmax <- X$q95
        xmean <- X$q50
        xlow <- X$q25
        xup  <- X$q75
        leg  <- c("q05-q95", "q25-q75")
      } else {
        xmin <- X$Min
        xmax <- X$Max
        xmean <- X$Mean
        xlow <- X$Mean-X$Sd
        xup  <- X$Mean+X$Sd
        leg <- c("Min-Max", "Mean+-sd")
      }
    
      # plotting parameters for deSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, ij)
      Dotpoints <- extractdots(dotpoints, ij)
      
      Xlog <- Ylog <- FALSE
      if (! is.null(Dotmain$log)) { 
        Ylog  <- length(grep("y",Dotmain$log))
        Xlog  <- length(grep("x",Dotmain$log))
      }       
      if (!xyswap) {
        if (is.null(yylim[[ij]]))
          Dotmain$ylim <- Range(NULL,c(xmin, xmax, obs[,io]), Ylog)
        else
          Dotmain$ylim <- yylim[[ij]]
        
        if (is.null(xxlim[[ij]]))
          Dotmain$xlim <- Range(NULL,X$x,Xlog)
        else
          Dotmain$xlim <- xxlim[[ij]]
        do.call("plot", c(alist(X$x, xmean, type = "n"), Dotmain))
        polygon(c(X$x, rev(X$x)), c(xmin, rev(xmax)), col = col[1], border = NA)
        polygon(c(X$x, rev(X$x)), c(xlow, rev(xup)), col = col[2], border = NA)
        do.call("lines", c(alist(X$x, xmean), Dotpoints) )
        if (! is.na(io))   
           for (j in 1: nobs) 
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0) 
                do.call("points", c(alist(obs[i.obs, 1], obs[i.obs, io]), 
                         extractdots(Obspar, j) ))        
        
      } else {
        if (is.null(xxlim[[ij]]))
          Dotmain$xlim <- range(c(xmin, xmax, obs[,io]), na.rm = TRUE)
        else
          Dotmain$xlim <- xxlim[[ij]]
  
        if (is.null(yylim[[ij]]))
          Dotmain$ylim <- rev(range(X$x))
        else
          Dotmain$ylim <- yylim[[ij]]
  
        do.call("plot", c(alist(xmean, X$x, type = "n"), Dotmain))
        polygon(c(xmin, rev(xmax)), c(X$x, rev(X$x)), col = col[1], border = NA)
        polygon(c(xlow, rev(xup)), c(X$x, rev(X$x)), col = col[2], border = NA)
        do.call("lines", c(alist(xmean, X$x), Dotpoints) )
        if (! is.na(io))   
           for (j in 1: nobs) 
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0) 
                do.call("points", c(alist(obs[i.obs, io], obs[i.obs, 1]), 
                         extractdots(Obspar, j) ))        
       }
    } # end i loop
    if (! is.null(legpos))
       legend(legpos, fill = c(col[1], col[2]),
               legend = leg, bty = "n")
  } else {            # nx =1 one summary per variable
    X <- x[Select, ]
    dotchart(X$Mean, labels = rownames(X), xlim = range(c(X$Min, X$Max)), ...)
    ## add ranges
    nr <- nrow(X)

    for (i in 1:nr) {
      segments(X$Min[i], i, X$Max[i], i, lty = 1)
      segments(X$q05[i], i, X$q95[i], i, lty = 1)
      segments(X$q25[i], i, X$q75[i], i, lwd = 3)
    }
    if (! is.null(legpos))
    legend(legpos, lwd = 1:3, legend = c("Min-Max", "q25-q75", "q05-q95"), bty = "n")
  }
}
