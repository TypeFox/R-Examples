plotGapfillCube <- function(
  ##title<< Plot an overview of a the results of a SSA gapfilling (from a ncdf file).
  file.orig                ##<< object to plot: file name or file.con object linking to the original (unfilled) ncdf file
  , file.filled = sub('[.]nc', '_gapfill.nc', file.orig)           ##<< object to plot: file name or file.con object linking to the filled ncdf file
  , file.prefill = ''      ##<< object to plot: file name or file.con object linking to the prefilled ncdf file                      
  , data.orig = c()        ##<< array: data (see file.orig). Can be supplied to prevent
                           ##          the reloading of huge datacubes.
  , data.filled = c()      ##<< see data.orig
  , data.prefill = c()     ##<< see data.orig                
  , n.series = 16          ##<< integer: how many example series to plot
  , lwd = 2                ##<< graphical parameter, see ?par
  , max.cores = 1          ##<< integer: amount of cores to use for parallelized computations.
  , ...
  )
  ##description<<
  ## This function plots some overview statistics of the results of a gapfilling
  ## (i.e. the results of a call to gapfillNcdf) .
  ##\if{html}{\out{<img src="../doc/visualize_ncdf_demo.png" alt="image ..visualize_ncdf_demo should be here"/>}}\ifelse{latex}{}{}
{
  ##TODO facilitate datacube input
  ##TODO include plotNLines capabilities
  
  ## preparation
  Sys.setenv(R_PARALLEL_PORT =  sample(49152:65535, 1))
  cl <- makeCluster(min(c(detectCores(), max.cores)))    
  con.orig   <- open.nc(file.orig)
  con.filled <- open.nc(file.filled)
  var.filled <- readNcdfVarName(file.filled)
  var.orig   <- readNcdfVarName(file.orig)

  ## load data
  if (length(data.orig) == 0) {
    cat('Loading original data ...\n')
    data.orig   <- transNcdfRotate(data.object = con.orig, var.name = var.orig)$data
  }
  if (length(data.filled) == 0) {
    cat('Loading gapfilled data ...\n')    
    data.filled <-  transNcdfRotate(data.object = con.filled, var.name = var.filled)$data 
  }
  if (length(data.prefill) == 0 & nchar(file.prefill)!=0) {
    cat('Loading pre-gapfilled data ...\n')
    con.prefill  <- open.nc(file.prefill)
    var.prefill  <- readNcdfVarName(file.prefill)
    data.prefill <-  transNcdfRotate(data.object = con.prefill, var.name = var.prefill)$data
  }

  
  
  ## calculate datacube info
  cat('Doing calculations ...\n')
  cube.info.orig     <- parApply(cl = cl, data.orig, 1:2, getVecInfo)       
  cube.info.filled   <- parApply(cl = cl, data.filled, 1:2, getVecInfo)
  stopCluster(cl)
  dimnames(cube.info.orig)[[1]] <- c('min', 'max', 'mean', 'sdev', 'range', 'ratio na',  'ratio na inner', 'ratio inf')
  dimnames(cube.info.filled)[[1]] <- c('min', 'max', 'mean', 'sdev', 'range', 'ratio na',  'ratio na inner', 'ratio inf')
  cube.info.agg       <- array(NA, dim=c(2, 9))
  dimnames(cube.info.agg) <- list(c('orig', 'filled'), c('min', 'max', 'mean', 'sd', 'range', 'ratio_full', 'ratio_empty', 'ratio_partial', 'ratio_continous'))
  for (dataset in c('orig', 'filled')) {
    data.t = get(paste('cube.info.', dataset, sep = ''))
    for (measure in c('min', 'max')) 
      cube.info.agg[dataset, measure] <- do.call(measure, list(data.t[measure, ,], na.rm = TRUE))
    for (measure in c('mean', 'sd')) 
      cube.info.agg[dataset, measure] <- do.call(measure, list(get(paste('data.', dataset, sep = '')), na.rm = TRUE))
    cube.info.agg[dataset, 'range']   <-   do.call('max', list(data.t['range', ,], na.rm = TRUE))
    cube.info.agg[dataset, 'ratio_empty']      <- sum(data.t['ratio na', , ]==1) / prod(dim(data.t['ratio na', , ]))
    cube.info.agg[dataset, 'ratio_full']       <- sum(data.t['ratio na', , ]==0) / prod(dim(data.t['ratio na', , ]))
    cube.info.agg[dataset, 'ratio_partial']    <- sum((data.t['ratio na', , ] > 0 & data.t['ratio na', , ] < 1)) / prod(dim(data.t['ratio na', , ]))
    cube.info.agg[dataset, 'ratio_continous']  <- sum(data.t['ratio na inner', , ]==0, na.rm = TRUE) / prod(dim(data.t['ratio na inner', , ]))
  }


  ## check consistency
  range.filled <- range(c(cube.info.filled['min', , ], cube.info.filled['max', , ]), na.rm = TRUE)
  range.orig   <- range(c(cube.info.orig['min', , ], cube.info.orig['max', , ]), na.rm = TRUE)

  if (sum((c(-1, 1) * range.orig) <= (c(-1 , 1) * range.filled)) != 2)   # if range orig exceeds range filled
    stop('Range orig is bigger than filled! Check code and input files!')
  if(sum(cube.info.orig['ratio na',,] - cube.info.filled['ratio na',,] < 0) > 0)
    stop('Filled dataset seems to have more missings than orig for some grids. Check code and input data!')
 
  
  ## do plots
  cat('Doing plots ...\n')
  grids.valid   <- which(cube.info.orig['ratio na', , ] < 1, arr.ind = TRUE)
  ind.rand      <- round(runif(16, 1, dim(grids.valid)[1]), digits = 0)
  ind.lat.orig  <- grids.valid[ind.rand, 1]
  ind.long.orig <- grids.valid[ind.rand, 2]
  ind.orig      <- matrix(NA, ncol = 3, nrow = length(ind.rand))
  colnames(ind.orig) <- c('lat', 'long', 'time')
  
  ind.orig[,'lat']  <- ind.lat.orig
  ind.orig[,'long'] <- ind.long.orig
  ind.orig[,'time']      <- '..'
  
  
  ## define color specs
  col.palette <- colorRampPalette(c('blue', 'yellow', 'red'))
  cols  <- rep(brewer.pal(ceiling(n.series/2), 'Set1'), times = 2)
  lty  <- rep(c(1), each = ceiling(n.series/2)*2)
  if (n.series > 8) {
    cols <- rep(cols, each = 2)
    lty  <- rep(c(1, 2), each = ceiling(n.series/2))
  }
  
  ## plot maps
  if(names(dev.cur()) == 'X11')
    dev.new()
  layout(matrix(c(1:9),byrow=TRUE,ncol=3),
         heights=c(1,1,1,1))
  par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(2, 0, 0, 2), oma = c(0, 2, 4, 0))
  pars.plot =  c('ratio na', 'min', 'max')
  for (i in 1:length(pars.plot)) {
    par.t = pars.plot[i]
    for (j in 1:3){
      if (j == 1) {
        data.t = cube.info.orig[par.t, , ]
      } else if (j == 2) {
        data.t = cube.info.filled[par.t, , ]
      } else if (j == 3) {
        data.t = cube.info.orig[par.t, , ] - cube.info.filled[par.t, , ]
      }
      if (j == 1)
        zlim <- range(c(cube.info.orig[par.t, , ], cube.info.filled[par.t, , ]), na.rm = TRUE)
      if (j == 3)
        zlim <- rangeZeroEqui(range(data.t, na.rm = TRUE))
      if (sum(!is.na(data.t)) > 0 ) {
        plotImageRotated(data.t,  xlab = '', xaxt = 'n', yaxt = 'n', col = col.palette(60),
                      zlim = range(pretty(zlim), scale = FALSE), useRaster = TRUE)
      } else {
        plot.new()
      }   
    }
  }
  labelMargins(pars.plot, las = 3, side = 2, outer = TRUE, cex= 2, line = .2)
  labelMargins(c('orig', 'filled', 'orig - filled'), side = 3, outer = TRUE, cex = 2, line =0.5)
  
  if(names(dev.cur()) == 'X11')
    dev.new()
  layout(matrix(c(1:2),byrow=TRUE,ncol=1),
         heights=c(1,1))
  par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(2, 0, 0, 2), oma = c(0, 2, 4, 0), xpd = FALSE)
  breaks = seq(min(cube.info.filled['min', , ], na.rm = TRUE),
               max(cube.info.filled['max', , ], na.rm = TRUE), length.out = 200)
  hst.orig    <- hist(data.orig, breaks = breaks, plot = FALSE)
  hst.filled  <- hist(data.filled, breaks = breaks, plot = FALSE)
  plotBG(rgb(.5,.5,.5))
  plot(hst.filled, xlim = range(c(hst.orig$mids, hst.filled$mids)),
       ylim = c(0, max( range(c(hst.orig$counts, hst.filled$counts)))),
       col = 'red', xlab = '', main = 'data')
  par(new=TRUE)  
  plot(hst.orig, xlim = range(c(hst.orig$mids, hst.filled$mids)),
       ylim = c(0, max( range(c(hst.orig$counts, hst.filled$counts)))), xlab = '',
       main = '', col = 'black')
  box()
  text(userCoords(c(0.8,0.9),c(0.9, 0.9)), labels =  c('filled', 'orig'),
       col = c('red', 'black'), cex = 2)
  logDensities <- log(c(hist(cube.info.filled['ratio na',,], breaks =seq(0,1,length.out=100), plot = FALSE)$density,
           hist(cube.info.orig['ratio na',,], breaks =seq(0,1,length.out=100), plot = FALSE)$density))
  yRange <- range(logDensities[is.finite(logDensities)])
  hst.filled  <- logHist(cube.info.filled['ratio na',,], breaks =100, col = 'red',
                         pch = 16, xlab = '', main = '', xlim = c(0,1), ylim = yRange)
  par(new = TRUE)
  hst.orig  <- logHist(cube.info.orig['ratio na',,], breaks =100, ylim = yRange,
          cex = 1.1, xlab = 'ratio of missing values per grid point', main = '', xlim = c(0,1))
  text(userCoords(c(0.7,0.9),c(0.9, 0.9)), labels =  c('filled', 'orig'),
       col = c('red', 'black'), cex = 2)
  mtext(side = 2 , text = 'log(density)', las = 3, line = 1)
   
  ## plot example series
  chars.plot <- c('ratio na', 'range', 'sdev', 'mean')
  if ( (length(data.prefill) != 0)) {
    chars.plot <- c(chars.plot, 'prefilling')
    ind.prefilled  <- is.na(data.prefill) & !is.na(data.orig)
  }
  for (characteristic in chars.plot) {
    if (characteristic == 'ratio na') {
      ratios.take <- seq(0, max(cube.info.orig['ratio na', , ][cube.info.filled['ratio na', , ] == min(cube.info.filled['ratio na', , ], na.rm = TRUE)], na.rm = TRUE), length.out = n.series)
      ind.valid   <- which(cube.info.filled['ratio na inner',,] == 0)
      ind.plot    <- indexVec2Matrix(ind.valid[whichClosest(ratios.take, cube.info.orig[characteristic, , ][ind.valid])], dim = dim(cube.info.orig)[-1])
    } else if (characteristic == 'prefilling') {
      ind.use        <- order(as.vector(apply(ind.prefilled, c(1,2), sum))/dim(ind.prefilled)[3], decreasing = TRUE)[1:n.series]      
      ind.plot    <- indexVec2Matrix(ind.use, dim = dim(cube.info.orig)[-1])
    } else {
      ind.valid    <- which(cube.info.filled['ratio na inner',,] == 0 & cube.info.filled['ratio na',,] < 0.95)
      ind.sorted   <- ind.valid[order(cube.info.filled[characteristic, , ][ind.valid], decreasing = TRUE)]
      if (is.element(characteristic, c( 'range', 'sdev')))
        ind.sorted <- ind.sorted[cube.info.filled[characteristic, , ][ind.sorted] !=0]
      ind.plot     <- indexVec2Matrix(ind.sorted[c(1:floor(n.series/2), (length(ind.sorted) - floor(n.series/2)):length(ind.sorted) )], dim = dim(cube.info.orig)[-1] )
    }
    if(names(dev.cur()) == 'X11')
      dev.new()
    layout(matrix(1:n.series, n.series, 1))
    par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2))
    args.extract =  c(list(data.filled), list(TRUE, TRUE, TRUE))
    args.orig    =  c(list(data.orig), list(TRUE, TRUE, TRUE))
    for (i in 1:n.series) {
      plotBG(rgb(0.9,0.9,0.9))
      args.extract[[2]] <- ind.plot[i,1]
      args.extract[[3]] <- ind.plot[i,2]
      args.orig[[2]]    <- ind.plot[i,1]
      args.orig[[3]]    <- ind.plot[i,2]
      y.data            <- do.call('[', args.extract)
      y.data.orig       <- do.call('[', args.orig)
      if (length(data.prefill) != 0) {
        ind.prefill  <- ind.prefilled[ind.plot[i,1], ind.plot[i,2], ]
      } else {
        ind.prefill     <- rep(FALSE, dim(data.filled)[3])
      }
      if (sum(!is.na(y.data)) > 0 ) {
        plot(y.data, col = 'red', pch = 16, type = 'b', lty = 2)
        text(userCoords(0.01,0.5), pos = 4, paste(unlist(args.extract[c(2,3)]), collapse = ','), cex = 2 )
        points(y.data.orig, col = 'black', pch = 16)
        if (sum(ind.prefill) > 0)
           points(which(ind.prefill), y.data.orig[ind.prefill], col = 'green', pch = 16)

      } else {
        plot.new()      
      }
    }
    text(userCoords(c(0.7,0.8, 0.9),c(0.4, 0.4, 0.4)), labels =  c('filled', 'orig', 'prefilled'),
        col = c('red', 'black', 'green'), cex = 2)
    mtext(characteristic, outer = TRUE, side = 3, cex = 2)
  }
 ## return stuff
  ##value<<
  ## some overview statistics of the different datacubes.
  invisible(cube.info.agg)
}

