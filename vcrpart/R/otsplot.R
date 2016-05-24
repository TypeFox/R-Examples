##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2014-11-10
##'
##' Description:
##' Ordinal time series plot and utility functions.
##'
##' Contents:
##'
##' otsplot_control:      control plot parameters
##' otsplot_filter:       creates filtering function for
##'                       fading out patterns
##' otsplot.default:      default function
##' print.otsplot:        print method for 'otsplot'
##' otsplot_panel:        panel function for otsplot
##' otsplot_data2str:     create trajectory strings
##' otsplot_layoutObjFun: layout optimization function
##' otsplot_layout:       layout optimizer
##' otsplot_color:        function to colour lines
##' linear:               linear colour filter function
##' minfreq:              minimal frequency for lines to filter
##' cumfreq:              filter a given proportion of most frequent
##'                       sequences
##'
##' Modifications:
##' 2014-11-10: added 'seed' argument to 'otsplot'
##' 2014-09-08: partial substitution of 'rep' by 'rep.int'
##' 2014-09-07: added header
##' -------------------------------------------------------- #

otsplot_control <- function(cex = 1, lwd = 1/4, col = NULL,
                            hide.col = grey(0.8), seed = NULL,
                            lorder = c("background", "foreground") ,
                            lcourse = c("upwards", "downwards"),
                            grid.scale = 1/5, grid.lwd = 1/2,
                            grid.fill =  grey(0.95), grid.col = grey(0.6),
                            layout = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
                            strip.fontsize = 12, strip.fill =  grey(0.9),
                            pop = TRUE, newpage = TRUE, maxit = 500L) {

  ## match arguments with duplicated entries in defaults
  lorder <- match.arg(lorder)
  lcourse <- match.arg(lcourse)

  ## checks
  stopifnot(is.logical(pop))
  stopifnot(is.logical(newpage))
  stopifnot(is.null(layout) | (is.numeric(layout) && length(layout) == 2L))
  stopifnot(is.numeric(margins) && length(margins) == 4L)
  stopifnot(is.numeric(maxit) && maxit > 0L)

  ## set list for plot parameters
  plot_gp <- list(cex = cex, lwd = lwd, col = col,
                  hide = hide.col, noobs = "black",
                  lorder = lorder, lcourse = lcourse,
                  alpha = 1, sf_cex = 0.5, sf_cex_leaves = 1)
  
  ## set list for translation zone plot parameters
  grid_gp <- list(scale = grid.scale, lwd = grid.lwd,
                  fill = grid.fill, col = grid.col)
  
  ## set list for strips plot parameters
  strip_gp <- gpar(fontsize = strip.fontsize,
                   fill = strip.fill)
   
  ## return a 'otsplot_control' object
  return(structure(
           list(plot_gp = plot_gp,
                grid_gp = grid_gp,
                strip_gp = strip_gp,
                layout = layout,
                margins = margins,
                pop = pop,
                newpage = newpage,
                maxit = maxit),
           class = "otsplot_control"))
}


otsplot_filter <- function(method = c("minfreq", "cumfreq", "linear"),
                           level = NULL) {
    method <- match.arg(method)
    if (is.null(level) && method %in% c("minfreq", "cumfreq"))
        stop("'otsplot_filter' requires an inpute for 'level'.")
    method <- switch(method,
                     linear = linear,
                     minfreq = minfreq,
                     cumfreq = cumfreq)
    return(structure(list(method = method, level = level),
                     class = "otsplot_filter"))
}


otsplot.default <- function(x, y, subject, weights, groups,
                            control = otsplot_control(), filter = NULL, 
                            main, xlab, ylab, xlim, ylim, ...) {  

  mc <- match.call()
  
  ## check filter argument
  if (is.character(filter))
    filter <- otsplot_filter(method = filter[1L], level = list(...)$level)
  if (!is.null(filter))
    stopifnot(inherits(filter, "otsplot_filter"))
  stopifnot(inherits(control, "otsplot_control"))
  
  ## check control argument
  stopifnot(inherits(control, "otsplot_control"))
  cArgs <- list(...)
  cArgs <- cArgs[intersect(names(cArgs), names(formals(otsplot_control)))]
  cArgs <- do.call("olmm_control", cArgs)
  control <- appendDefArgs(cArgs, control)

  ## set seed
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  if (!is.null(control$seed)) set.seed(control$seed)
  RNGstate <- .Random.seed
  
  ## check compatibility between 'filter' and 'lorder' control argument
  if (!is.null(filter)) control$plot_gp$lorder <- "foreground"
  
  ## check and prepare raw data

  if (missing(weights)) weights <- rep.int(1.0, length(x))
  if (missing(groups)) groups <- factor(rep.int(1L, length(x)))
  if (!all(c(length(y), length(subject), length(groups), length(weights))
           == length(x)))
    stop("input vectors have different lengths.")
  
  ## check the x variable
  if (is.numeric(x)) {
    if (length(unique(x)) > 1L) {
      xdiff <- diff(sort(unique(x)))
      if (sum((xdiff / min(xdiff) -
               round(xdiff / min(xdiff),0)) != 0) == 0) {
        x <- factor(x, levels = seq(min(x, na.rm = TRUE),
                         max(x, na.rm = TRUE), min(xdiff)))
      } else {
        warning("problems with distances between 'x' positions. ",
                "The 'x' positions will not be illustrated adequately.")
        x <- factor(x)
      }
    } else { x <- factor(x) }  # pseudo case
  }
  
  ## y must be an ordered factor
  if (!is.factor(y)) stop("'y' must be a 'factor'.") 
  
  ## convert subject variable to factor
  if (!is.factor(subject)) subject <- factor(subject)

  ## remove NAs and omit.levels in x, y or subject
  subs <- is.na(subject) | is.na(x) | is.na(y) | is.na(weights)
  if (sum(subs) > 0L) {
    subject <- subject[!subs]
    x <- x[!subs]
    y <- y[!subs]
    weights <- weights[!subs]
    groups <- groups[!subs]
  }
  subject <- droplevels(subject)
  y <- droplevels(y)
  x <- droplevels(x)
  
  ## construct/ convert groups variable
  if (!is.factor(groups)) groups <- factor(groups)
  groups <- droplevels(groups)
  if (sum(is.na(groups)) > 1L) {
    groups <- factor(groups, levels = c(levels(groups), "unknown group"))
    groups[is.na(groups)] <- "unknown group"
  }
  groups <- unique(data.frame(subject, groups))[, 2]    
  if (length(groups) != nlevels(subject)) 
    stop("cannot link 'groups' with 'subject'.")        
  names(groups) <- levels(subject)
  
  ## construct/ convert weights variable
  if (!is.numeric(weights)) stop("'weights' are not numeric.")
  if (min(weights, na.rm = TRUE) < 0L) stop("negative 'weights'.")
  if (length(weights) == length(subject))
    weights <- unique(data.frame(subject, weights))[, 2L]
  if (length(weights) != nlevels(subject))
    stop(paste("cannot link 'weights' with 'subject'."))
  
  ## standardize weights
  wsubjectGroups.unscaled <- tapply(weights, groups, sum)
  tmp1 <- table(groups)[as.integer(groups)]
  tmp2 <- tapply(weights, groups, sum)[as.integer(groups)]
  weights <- weights / tmp2 * tmp1
  names(weights) <- levels(subject)
  
  ## subject variable
  nSubject <- nlevels(subject)
  subject <- as.integer(subject)
  
  ## x variable
  xLevs <- levels(x)
  nXLevs <- nlevels(x)
  x <- as.integer(x)
  
  ## y variable
  yLevs <- levels(y)
  nYLevs <- nlevels(y)
  y <- as.integer(y)
  
  ## groups variable
  nSubjectGroups <- table(groups)
  wsubjectGroups <- tapply(weights, groups, sum)
  groupsLevs <- levels(groups)
  nGroups <- nlevels(groups)
  groups <- as.integer(groups)
    
  ## order the data
  subs <- order(subject, x, y)
  subject <- subject[subs]
  x <- x[subs]
  y <- y[subs]
  
  ## Find curves to plot
    
  ## some needed variables
  x.list <- tapply(X=x,
                   INDEX=list(subject),
                   FUN=function(x){return(x)})  
  trajString <- tapply(X = y, INDEX = list(subject, x), paste,
                        collapse = ",")
  trajString <- apply(X = trajString, MARGIN = 1L, FUN = otsplot_data2str)
  
  ## find trajectories
  trajGroup <- as.integer(factor(trajString))
  trajSubject <- 1L:nSubject
  groups.name <- 1L:max(trajGroup)
  groups.subject <- tapply(trajSubject, trajGroup, function(x){x[1L]}) 
  
  ## prepare point data
    
  ## setup pts object
  pts <- data.frame(subject = subject[subject %in% groups.subject],
                    x = x[subject %in% groups.subject],
                    y = y[subject %in% groups.subject])
  pts$traj <- as.integer(factor(pts$subject,
                                levels = groups.subject,
                                labels = groups.name))
  ntraj <- max(pts$traj)
  
  ## if lines should to run from the highest to lowest y category
  ## in case of simultaneous observations
  if (control$plot_gp$lcourse == "downwards") {
    pts <- pts[order(pts$traj, pts$x,
                     factor(pts$y, levels = rev(sort(unique(pts$y))))), ,
               drop = FALSE]
  } 
  
  ## find multiple equal simultaneous observations
  ## and note their frequency
  tmp <- as.data.frame(table(x=pts$x, y=pts$y, traj = pts$traj),
                       responseName = "duplicates")    
  tmp <- tmp[tmp$duplicates != 0L,]
  pts <- merge(x = unique(pts), y = tmp, by = c("x", "y", "traj"),
               sort = FALSE)
  
  ## some useful column extractors
  freqCols <- paste("n", 1L:nGroups, sep = ".")
  widthCols <- paste("wdt", 1L:nGroups, sep = ".")
  colCols <- paste("col", 1L:nGroups, sep = ".")
   
  ## weighted point sizes
  pts[, freqCols] <- tapply(weights, list(trajGroup, groups), sum)[pts$traj, ]
  pts[is.na(pts)] <- 0.0
  
  ## jitter coordinates
    
  ## needed variables
  trajWeightsMax <- apply(X = pts[,freqCols, drop = FALSE], MARGIN = 2L,
                          FUN = tapply, list(pts$traj), max)
  trajProp <- scale(x = matrix(trajWeightsMax, ncol = nGroups),
                       center = FALSE, scale = wsubjectGroups)
  trajPropMax <- apply(X = trajProp, MARGIN = 1L, FUN = max)
  trajord <- order(trajPropMax, decreasing = TRUE)
  
  otsplot_jitter <- function(data, maxit) {
    
    ## create initial grid
    ngrid <- ngrid0 <- 10
    sl <- ceiling(ngrid * sqrt(trajPropMax)) # size
    ngrid <- ceiling(sqrt(sum(sl ^ 2)))
    grid <- matrix(0, ngrid, ngrid) # generate initial grid
    gridsubs <- 1L:(ngrid ^ 2) # identifiers
    data$xpos <- data$ypos <- NA      
    potpos <- c()
    blacklist <- c()
    
    ## the algorithm
    for (i in 1L:length(trajord)) {
      subscripts <- data$traj == trajord[i]
      count <- 0
      found <- FALSE
      while ((!found)&(count <= maxit)) {
        if ((i > 1L) & (count == 0)) {
          if (sl[trajord[i]] < sl[trajord[i-1L]]) {
            potpos <- c()
            blacklist <- c()
          }
        }
        count2 <- 0
        while (length(potpos) == 0) { # no potential position available
          if (count2 > 0) { # enlarge grid
            ngrid <- ngrid + 1L
            grid <- rbind(grid, rep.int(0, ngrid - 1L))
            grid <- cbind(grid, rep.int(0, ngrid))
            gridsubs <- seq(1L, ngrid^2, 1L)
            ## correct blacklist
            blacklist <-  blacklist + floor(blacklist/ngrid)
          }
          ## determine potential grid positions
          potpos <- (grid[gridsubs] == 0) & # occupied positions
          ((gridsubs %% ngrid) < (ngrid - sl[trajord[i]] + 2)) & # to close top
          (gridsubs < (ngrid * (ngrid - sl[trajord[i]] + 1L))) # to close left
          if (sl[trajord[i]] > 1L) { # filter for groups with ns > 1 field
            potpos <- potpos & (gridsubs %% ngrid != 0)
          }
          potpos <- which(potpos)
          potpos <- potpos[!potpos %in% blacklist]
          count2 <- count2+1L
        }     
        xpotpos <- floor(potpos / ngrid)+1L
        ypotpos <- potpos %% ngrid
        ypotpos[ypotpos == 0] <- ngrid
        posind <- sample(1L:length(potpos), 1L) # random assignement
        pos <- potpos[posind]
        xpos <- xpotpos[posind]
        ypos <- ypotpos[posind]
        xsubs <- seq(xpos, xpos + sl[trajord[i]] - 1L, 1L)
        ysubs <- seq(ypos, ypos + sl[trajord[i]] - 1L, 1L)
        if (sum(c(grid[ysubs, xsubs]) != 0) == 0) { # assign position
          data$xpos[subscripts] <- xpos
          data$ypos[subscripts] <- ypos
          grid[ysubs,xsubs] <- trajord[i] # register assignement in matrix
          potpos <- potpos[!potpos %in% which(grid == trajord[i])] 
          found <- TRUE
        } else { # nothing found
          count <- count + 1L
          blacklist <- c(blacklist, pos)
          potpos <- potpos[-posind]
        }
        if (count == maxit) {
          tmp1 <- data[data$traj == i,]
          tmp2 <- otsplot_data2str(y = tmp$y1, x = tmp$x1)
          warning(paste(" found no grid position for trajectory: ",
                        tmp2 ,", omit!",sep=""))
        }
      }
    }
    data$xjitter <- (data$xpos - 1L) / ngrid
    data$yjitter <- (data$ypos - 1L) / ngrid
    ## determine central plot coordinates
    
    return(list(jitter = data[,c("xjitter", "yjitter")],
                ngrid0 = ngrid0, ngrid = ngrid))
  }
  tmp <- otsplot_jitter(pts, control$maxit)
  pts[,c("xjitter","yjitter")] <- tmp$jitter
  ngrid0 <- tmp$ngrid0
  ngrid <- tmp$ngrid
  
  ## curve colors
  if (is.null(control$plot_gp$col)) { # colors are not specified by the user
    control$plot_gp$col <- rainbowPalette[rep.int(seq(1L, 8), ceiling(ntraj / 8))][1L:ntraj]
  } else {
    control$plot_gp$col <- rep_len(control$plot_gp$col, ntraj) 
  }
  control$plot_gp$col <- control$plot_gp$col[order(trajord)]
  if (is.character(control$plot_gp$col)) {
    control$plot_gp$col <- col2rgb(control$plot_gp$col)
    control$plot_gp$col <- rgb(control$plot_gp$col[1L,],
                               control$plot_gp$col[2L,],
                               control$plot_gp$col[3L,],
                               control$plot_gp$alpha * 255, maxColorValue = 255)
  }
  
  pts[, colCols] <- matrix(rep.int(control$plot_gp$col[pts$traj], nGroups), ncol = nGroups)
  control$plot_gp$noobs <- rep.int(control$plot_gp$noobs, nGroups)
  
  ## color gradients
  if (!is.null(filter)) {
        
    tmp1 <- matrix(NA, ncol = nGroups, nrow = ntraj)
    for (i in 1L:ncol(trajProp))
      tmp1[, i] <- filter$method(trajProp[,i], level = filter$level)
    tmp3 <- tmp1
    for (i in 1L:nrow(tmp1)) {
      for (j in 1L:ncol(tmp1)) {
        tmp3[i, j] <- otsplot_color(c(0, tmp1[i, j], max(tmp1[, j])), control$plot_gp$hide, c(control$plot_gp$col, control$plot_gp$noobs[1L])[i])[2]
      }
    }
    for (i in 1L:nGroups) {
      pts[, colCols[i]] <- tmp3[-nrow(tmp3), i][pts$traj]
    }
    control$plot_gp$noobs <- tmp3[nrow(tmp3), ]
  }
  
  ## xlab and ylab
  if (missing(main)) main <- NULL
  if (missing(xlab)) xlab <- deparse(match.call()$x)
  if (missing(ylab)) ylab <- NULL
  
  ## xlim and ylims
  if (missing(xlim)) {
    xlim <- c(1L - sqrt(control$grid_gp$scale) / 2,
              nXLevs + sqrt(control$grid_gp$scale) / 2)
    xlim <- xlim + c(-1, 1) * 0.025 * diff(range(xlim))
  }
  if (missing(ylim)) {
    ylim <- c(1 - sqrt(control$grid_gp$scale) / 2,
              nYLevs + sqrt(control$grid_gp$scale) / 2)
    ylim <- ylim + c(-1, 1) * 0.025 * diff(range(ylim))
  }
  
  ## point widths
  pts[, widthCols] <- data.frame(sqrt(scale(x = pts[, freqCols, drop = FALSE], center = FALSE, scale = wsubjectGroups)) * sqrt(control$grid_gp$scale) * ngrid0 / ngrid)
  
  ## data frame for line segments
  ntraj <- max(pts$traj)
  lns <- NULL
  for (i in unique(pts$traj)) {
    subs <- which(pts$traj == i)
    nsubs <- length(subs)
    if (length(subs) > 1L) {
      tmp <- data.frame(traj = rep.int(i, nsubs-1L),
                        x0 = pts$x[subs[-nsubs]],
                        y0 = pts$y[subs[-nsubs]],
                        x1 = pts$x[subs[-1L]],
                        y1 = pts$y[subs[-1L]])
      lns <- rbind(lns,tmp)
    }
  }
  
  if (!is.null(lns)) {
  
    ## merge coordinates and colors
    lns <- merge(x = lns, y = pts[, c("traj", "x", "y", "xjitter", "yjitter")], by.x = c("traj", "x0", "y0"), by.y = c("traj", "x", "y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
    lns <- merge(x = lns, y = pts[,c("traj","x","y","xjitter","yjitter", freqCols, colCols, widthCols)], by.x = c("traj","x1","y1"), by.y = c("traj","x","y"), all.x = TRUE, all.y = FALSE, sort = FALSE)
    colnames(lns)[6:9] <- c("x0jitter", "y0jitter", "x1jitter", "y1jitter")
  }
  
  if (nGroups > 1L) { # plot layout
    if (is.null(control$layout)) {
      optpr <- otsplot_layout(nGroups, nGroups, nGroups, 1L)
      nLayCols <- optpr[1L]
      nLayRows <- optpr[2L]
    } else {
      nLayCols <- control$layout[2L]
      nLayRows <- control$layout[1L]
    }
  } else { nLayCols <- 1L; nLayRows <- 1L }
  
  ## coordinates for background rectangles
  backrect <- expand.grid(xgrid = 1L:nXLevs, ygrid = 1L:nYLevs)

  pts <- reshape(pts, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1L:nGroups, direction = "long")
  lns <- reshape(lns, varying = list(freq = freqCols, width = widthCols, col = colCols), v.names = c("freq", "width", "col"), timevar = "groups", times = 1L:nGroups, direction = "long")

  ## reset seed
  assign(".Random.seed", oldSeed, envir=globalenv())
  
  ## plot data object
  return(structure(
           list(pts = pts, lns = lns, backrect = backrect, 
                trajProp = trajProp, trajPropMax = trajPropMax,
                ntraj = ntraj,
                nGroups = nGroups, groupsLevs = groupsLevs,
                nXLevs = nXLevs, nYLevs = nYLevs,
                xLevs = xLevs, yLevs = yLevs,
                main = main, xlab = xlab, ylab = ylab,
                xlim = xlim, ylim = ylim,
                nLayCols = nLayCols, nLayRows = nLayRows,
                ngrid = ngrid, ngrid0 = ngrid0,
                grid_gp = gpar(col = control$grid_gp$col,
                  lwd = control$grid_gp$lwd,
                  fill = control$grid_gp$fill),
                grid_scale = control$grid_gp$scale, 
                strip_gp = control$strip_gp,
                cex = control$plot_gp$cex,
                lwd = control$plot_gp$lwd, lorder = control$plot_gp$lorder,
                col_noobs = control$plot_gp$noobs,
                sf_cex = control$plot_gp$sf_cex,
                sf_cex_leaves = control$plot_gp$sf_cex_leaves,
                newpage = control$newpage, margins = control$margins,
                pop = control$pop), class = "otsplot"))
}


print.otsplot <- function(x, ...) {

  if (x$newpage) grid.newpage()

  axTicks <- pretty(1L:x$nXLevs)
  axTicks <- intersect(axTicks, 1L:x$nXLevs)
  
  if (x$nGroups == 1L) {
    
    pushViewport(plotViewport(xscale = x$xlim,
                              yscale = x$ylim,
                              default.units = "native",
                              margins = x$margins,
                              name = "otsplot"))
    grid.rect()
    otsplot_panel(x, 1L)
    grid.xaxis(axTicks, x$xLevs[axTicks])
    grid.yaxis(1L:x$nYLevs, x$yLevs)
    grid.text(x$xlab, y = unit(-3, "lines"))
    grid.text(x$ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"),
              gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
    
  } else {
  
    strUnit <- unit(2, "strheight", "A")
    pushViewport(plotViewport(layout = grid.layout(x$nLayRows, x$nLayCols),
                              name = x$name, margins = x$margins))
    
    for (i in 1L:x$nGroups) {

      ## cell position
      posr <- floor((i - 1L) / x$nLayCols) + 1L
      posc <- i - floor((i - 1L) / x$nLayCols) * x$nLayCols
      
      pushViewport(viewport(layout.pos.col = posc,
                            layout.pos.row = posr,
                            name = paste("otsplot_cell", i, sep = "_"),
                            layout = grid.layout(2, 1L,
                              heights = unit.c(strUnit, unit(1, "npc") - strUnit))))  
      
      ## header
      pushViewport(viewport(layout.pos.row = 1L, name = paste("strip", i , sep = ".")))
      grid.rect(gp = x$strip_gp)
      grid.text(x$groupsLevs[i], gp = x$strip_gp)
      upViewport(1L)

      ## plot
      pushViewport(viewport(layout.pos.row = 2L, name = paste("plot", i, sep = "_")))
      grid.rect()
      pushViewport(viewport(xscale = x$xlim,
                            yscale = x$ylim,
                            default.units = "native"))
      
      otsplot_panel(x, i)

      if (posc == 1L) grid.yaxis(1L:x$nYLevs, x$yLevs)
      if (i > x$nGroups - x$nLayCols) grid.xaxis(axTicks, x$xLevs[axTicks])
      
      upViewport(3L)
    }
    grid.text(x$xlab, y = unit(-3, "lines"))
    grid.text(x$ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(x$main, y = unit(1, "npc") + unit(2, "lines"),
              gp = gpar(fontface = "bold"))
    
    if (x$pop) popViewport(1L) else upViewport(1L)
  }
}


otsplot_panel <- function(x, groups) {

  xf <- diff(x$xlim) / as.numeric(convertWidth(unit(1, "npc"), "inches"))
  yf <- diff(x$ylim) / as.numeric(convertHeight(unit(1, "npc"), "inches"))
  tmp <- max(xf, yf)
  xf <- xf/tmp; yf <- yf/tmp
  
  ## extract data
  x$pts <- x$pts[x$pts$groups == groups, ]
  x$lns <- x$lns[x$lns$groups == groups, ]
    
  ## temporary adaption of line widths
  x$lns$lwd <- 96 * as.numeric(convertUnit(unit(x$lns$width, "native"), "inches", ifelse(xf < yf, "x", "y"), "dimension")) * min(c(xf, yf)) * x$cex * x$lwd
  
  ## translation zones
  grid.rect(x = unit(x$backrect$xgrid, "native"),
            y = unit(x$backrect$ygrid, "native"),
            width = unit(xf * sqrt(x$grid_scale), "native"),
            height = unit(yf * sqrt(x$grid_scale), "native"),
            gp = x$grid_gp)
  
  for (i in order(x$trajProp[,groups],
                  decreasing = (x$lorder == "background"))) {
    
    ## extract subscripts
    subspts <- (x$pts$traj %in% i) & (x$pts$freq > 0L)
    subssun <- subspts & (x$pts$duplicates > 1L)
    if (!is.null(x$lns)) {
      subslns <- (x$lns$traj %in% i) & (x$lns$freq > 0L)
      } else {
        subslns <- FALSE
      }    
    
      ## points
    if (sum(subspts) > 0L) {
      grid.rect(x = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[subspts], "native"), y = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[subspts], "native"), width = unit((x$pts$width * xf * x$cex)[subspts], "native"), height = unit((x$pts$width * yf * x$cex)[subspts], "native"), gp = gpar(col = x$pts$col[subspts], fill = x$pts$col[subspts], lwd = 0, lex = 0))
      }

    ## lines
    if (sum(subslns) > 0L) {
      grid.segments(x0 = unit((x$lns$x0 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[subslns], "native"), y0 = unit((x$lns$y0 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[subslns], "native"), x1 = unit((x$lns$x1 + xf * (- 0.5 + x$lns$x0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[subslns], "native"), y1 = unit((x$lns$y1 + yf * (- 0.5 + x$lns$y0jitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$lns$traj] / 2))[subslns], "native"), gp = gpar(lwd = x$lns$lwd[subslns], col = x$lns$col[subslns], lineend = "round"))
    }
    
    ## sunflowers
    if (sum(subssun) > 0L) {
      
      grid.points(x = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[subssun], "native"), y = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[subssun], "native"), gp = gpar(cex = x$sf_cex), pch = 16)
      i.multi <- which(subssun) # stolen from sunflowerplot()
      ppin <- par("pin")
      pusr <- par("usr")
      xr <- x$pts$width[subssun] / 2 * x$sf_cex_leaves * xf
      yr <- x$pts$width[subssun] / 2 * x$sf_cex_leaves * yf
      i.rep <- rep.int(i.multi, x$pts$duplicates[subssun])
      z <- numeric()
      for (k in i.multi) {
        z <- c(z, 1:x$pts$duplicates[k])
      }
      deg <- (2 * pi * z)/x$pts$duplicates[i.rep]
      
      grid.segments(x0 = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], "native"), y0 = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep], "native"), x1 = unit((x$pts$x + xf * (- 0.5 + x$pts$xjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + xr * sin(deg), "native"), y1 = unit((x$pts$y + yf * (- 0.5 + x$pts$yjitter * sqrt(x$grid_scale) + 1/2 * (1 - sqrt(x$grid_scale)) + c(sqrt(x$trajPropMax) * sqrt(x$grid_scale) * x$cex * x$ngrid0 / x$ngrid)[x$pts$traj] / 2))[i.rep] + yr * cos(deg), "native"))
    }
  }
}
  

otsplot_data2str <- function(y, x = 1L:length(y),
                         pa.left = "(", pa.right = ")",
                         pot= "^" ,con = "-") {
  subscripts <- !is.na(y);
  ret <- paste(pa.left, y, pa.right, pot, x, sep = "");
  ret <- ret[subscripts];
  ret <- paste(ret, collapse = con)
  return(ret)
}


otsplot_layoutObjFun <- function(nXLevs, nYLevs, n, c = 0.5) {
  control <- nXLevs * nYLevs - n # number of empty windows
  ## optimization value: function of number of empty windows
  ## and ncol/nrow ratio 
  opt <- control / n + c * abs(1 - min(nXLevs, nYLevs) / max(nXLevs, nYLevs))
  return(c(control,opt))
}


otsplot_layout <- function(nXLevs, nYLevs, npanels, c = 1L) {
  minvalue1 <- otsplot_layoutObjFun(nXLevs - 1L, nYLevs, npanels, c)
  minvalue2 <- otsplot_layoutObjFun(nXLevs, nYLevs - 1L, npanels, c)
  while (minvalue1[1L] >= 0 | minvalue2[1L] >= 0) {
    if ((minvalue1[1L] >= 0) & (minvalue2[1L] >= 0)) {
      if (minvalue1[2] < minvalue2[2]) {
        nXLevs <- nXLevs - 1L } else { nYLevs <- nYLevs - 1L }
    } else {
      if (minvalue1[2] >= 0) {
        nXLevs <- nXLevs - 1L
      } else {
        if (minvalue2[2] >= 0) {
          nYLevs <- nYLevs - 1L
        }
      }}
    minvalue1 <- otsplot_layoutObjFun(nXLevs - 1L, nYLevs, npanels, c)
    minvalue2 <- otsplot_layoutObjFun(nXLevs, nYLevs - 1L, npanels, c)
  }
  return(c(nXLevs, nYLevs))
}


otsplot_color <- function(value, col1, col2) {
  mp <- colorRamp(c(col1, col2))
  col <- rgb(mp(value), maxColorValue = 255)
}


linear <- function(x, ...) {
  return((x - min(x)) / diff(range(x)))
}


minfreq <- function(x, level = 0.05) {
  return(1*(x >= level))
}


cumfreq <- function(x, level = 0.75) {
  tmp <- which(cumsum(sort(x, decreasing = TRUE)) >= level)[1L]
  ret <- vector("logical", length(x))
  ret[order(x, decreasing = TRUE)[1L:tmp]] <- TRUE
  return(1 * ret)
}
