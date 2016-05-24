plotDatacube <- function(
  ##title<< Visualize/plot an overview of a netCDF file
  data.object         ##<< object to plot: file name or file.con object linking to a netCDF file
  , data = c()        ##<< array: data to plot. Can be passed to the function to
                      ##   prevent the repeated loading of huge netCDF data.
  , fourth.dim = 0     ##<< position in possible forth dimension (height, spectral band etc) to plot
  , var.name = 'auto' ##<< character string: name of the variable to plot
  , parallel = FALSE  ##<< logical: Whether to parallelize the computations.
  , max.cores = 16    ##<< integer: maximum amount of cores to use for the
                      ##   parallelized computations.
  , n.series = 16     ##<< integer: amount of example series to plot
  , lwd = 2           ##<< integer: graphical parameter, see ?par                     
  , ...
)
  ##description<<
  ## This function plots some overview statistics of a netCDF file.
  ##\if{html}{\out{<img src="../doc/visualize_ncdf_demo.png" alt="image ..visualize_ncdf_demo should be here"/>}}\ifelse{latex}{}{}
{
  ##TODO facilitade datacube input
  ##TODO include plotNLines capabilities

  cat('Preparing stuff ...\n')
  if (parallel) {
    cat('Registering cluster ...\n')
    cl <- makeCluster(max.cores)
  }  
  if (class(data.object) == 'NetCDF') {
    file.con   = data.object
  } else {
    file.con   = open.nc(data.object)
  }
  
  if (var.name == 'auto')
    var.name <- readNcdfVarName(file.con)
  
  lat.name   <- infoNcdfDims(file.con)$name[pmatch('lat', infoNcdfDims(file.con)$name)]
  long.name  <- infoNcdfDims(file.con)$name[pmatch('lon', infoNcdfDims(file.con)$name)]
  
  dim.long   <- match(long.name, infoNcdfDims(file.con)$name)
  dim.lat    <- match(lat.name, infoNcdfDims(file.con)$name)
  dim.time   <- match('time', infoNcdfDims(file.con)$name)
  if (is.element(long.name, infoNcdfVars(file.con,dimvars=TRUE)[, 'name'])) {
    longitudes <- round(var.get.nc(file.con, long.name), digits = 0)
  } else {
    longitudes <-  1:infoNcdfDims(file.con)[dim.long, 'length']
  }
  if (is.element(long.name, infoNcdfVars(file.con,dimvars=TRUE)[, 'name'])) {
    latitudes  <- round(var.get.nc(file.con, lat.name), digits = 0)  
  } else {
    latitudes  <-  1:infoNcdfDims(file.con)[dim.lat, 'length']
  }
  if (is.element('time', infoNcdfVars(file.con,dimvars=TRUE)[, 'name'])) {
    time       <- convertDateNcdf2R(file.con)
  } else {
    time  <-  1:infoNcdfDims(file.con)[infoNcdfDims(file.con)$name == 'time', 'length']
  }
  dims.par   <- infoNcdfDims(file.con)$name[var.inq.nc(file.con, var.name)$dimids + 1]
  
  ## sort dataframe
  aperm.data <- pmatch(c('lat', 'lon', 'time'), dims.par)
  if (length(data) == 0) {
    cat('Loading data ...\n')
    data       <- var.get.nc(file.con, var.name)
  }
  cat('Transposing datacube ...\n')
  
  if (length(fourth.dim) > 1 || fourth.dim != 0) {
    aperm.data <- c(aperm.data, 4)
    length.fourth.dim <- dim(data)[4]
  } else {
    length.fourth.dim <- 1
  }
  data.cube.unsort   <- aperm(data, aperm.data)
  
  if (length(fourth.dim) > 1 || fourth.dim != 0) {
    data.cube.sort   <- data.cube.unsort[order(latitudes, decreasing = TRUE), order(longitudes), ,]
  } else {
    data.cube.sort   <- data.cube.unsort[order(latitudes, decreasing = TRUE), order(longitudes),]
  }
  data.cube.sort     <- array(data.cube.sort, dim = c(dim(data.cube.sort)[1:3], length.fourth.dim) ) 
  
  ## calculate datacube info
  cat('Doing calculations ...\n')
  if (parallel) {
    cube.info          <- parApply(cl, data.cube.sort, c(1,2,4), getVecInfo)       
  } else {
    cube.info          <- apply(data.cube.sort, c(1,2,4), getVecInfo)
  } 
#  if  (length(fourth.dim) == 1 && fourth.dim == 0)
#    cube.info.t              <- array(cube.info, dim = c(dim(cube.info)[1:2], 1, dim(cube.info)[3]))

  if (dim(cube.info)[4] > 1) {
    cube.info.agg      <- data.frame(min = apply(cube.info['min',,,], 3, min, na.rm = TRUE),
                                     max = apply(cube.info['max',,,], 3, max, na.rm = TRUE),
                                     mean = apply(cube.info['mean',,,], 3, mean, na.rm = TRUE),
                                     ratio.na = apply(cube.info['ratio na',,,], 3, mean, na.rm = TRUE),
                                     series.empty = apply(cube.info['ratio na',,,], 3, function(x) sum(x == 1)) / prod(dim(data.cube.sort)[2:3]),
                                     series.full = apply(cube.info['ratio na',,,], 3, function(x) sum(x == 0)) / prod(dim(data.cube.sort)[2:3]), 
                                     series.gappy = apply(cube.info['ratio na',,,], 3, function(x) sum(x != 0 & x != 1)) / prod(dim(data.cube.sort)[2:3]))
  } else {
    cube.info.agg      <- data.frame(min = min(cube.info['min',,,1], na.rm = TRUE),
                                     max = max(cube.info['max',,,1], na.rm = TRUE),
                                     mean = mean(cube.info['mean',,,1], na.rm = TRUE),
                                     ratio.na = mean(cube.info['ratio na',,,1], na.rm = TRUE),
                                     series.empty = sum(cube.info['ratio na',,,] == 1 ) / prod(dim(data.cube.sort)[2:3]),
                                     series.full = sum(cube.info['ratio na',,,] == 0) / prod(dim(data.cube.sort)[2:3]), 
                                     series.gappy = sum(cube.info['ratio na',,,] != 0 & cube.info['ratio na',,,] != 1) / prod(dim(data.cube.sort)[2:3]))
  }
  cat('Doing plots ...\n')
  for (h in 1:length(fourth.dim)) {
    fourth.dim.t = fourth.dim[h]
    if (length(fourth.dim) == 1 && fourth.dim == 0)
      fourth.dim.t <- 1
    
    grids.valid <- which(cube.info['ratio na', , ,fourth.dim.t] < 1, arr.ind = TRUE)
    ind.rand      <- round(runif(16, 1, dim(grids.valid)[1]), digits = 0)
    ind.lat.orig  <- (1:length(latitudes))[order(latitudes, decreasing = TRUE)][grids.valid[ind.rand, 1]]
    ind.long.orig <- (1:length(longitudes))[order(longitudes)][grids.valid[ind.rand, 2]]
    ind.orig      <- matrix(NA, ncol = {if(length(fourth.dim) == 1 && fourth.dim == 0){3} else {4}}, nrow = length(ind.rand))
    if  (length(fourth.dim) == 1 && fourth.dim == 0) {
      colnames(ind.orig) <- c('lat', 'long', 'time')
    } else {
      colnames(ind.orig) <- c('lat', 'long', 'time', 'fourth.dim')
      ind.orig[,4]           <- fourth.dim.t
    }      
    ind.orig[,'lat']  <- ind.lat.orig
    ind.orig[,'long'] <- ind.long.orig
    ind.orig[,'time']      <- '..'
    if (length(fourth.dim) > 1 || fourth.dim != 0)
      ind.orig[,4]           <- fourth.dim.t
    
    ## define color specs
    col.palette <- colorRampPalette(c('blue', 'yellow', 'red'))
    cols  <- rep(brewer.pal(ceiling(n.series/2), 'Set1'), times = 2)
    lty  <- rep(c(1), each = ceiling(n.series/2)*2)
    if (n.series > 8) {
      cols <- rep(cols, each = 2)
      lty  <- rep(c(1, 2), each = ceiling(n.series/2))
    }
    
    ## plot maps
    z.range     <- range(cube.info[, , , fourth.dim.t],na.rm = TRUE, finite = TRUE)
    if (length(fourth.dim) > 1 & names(dev.cur()) == 'null device' )
      dev.new()        
    layout(matrix(c(2,4,6,8,10,12,1,3,5,7,9,11,13,13,13,13,13,13),byrow=TRUE,ncol=6),
        heights=c(0.1,1,1))
    par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(2, 0, 0, 1), oma = c(0, 2, 4, 1))
    pars.plot = c('min', 'max', 'mean', 'sdev', 'ratio na', 'ratio inf')
    for (i in 1:length(pars.plot)) {
      if (sum(!is.na(cube.info[pars.plot[i], , , fourth.dim.t])) > 0 ) {
        plotImageRotated(cube.info[pars.plot[i], , , fourth.dim.t], row.vals = latitudes,
            col.vals = longitudes, xlab = '', col = col.palette(60),
            zlim = range(pretty(cube.info[pars.plot[i], , , fourth.dim.t])), scale = FALSE)
        if (i == 1) {
          y <-  latitudes[order(latitudes, decreasing=TRUE)][grids.valid[ind.rand[1:n.series], 1]]
          x <-  longitudes[grids.valid[ind.rand[1:n.series], 2]]
          points(x, y, cex = 4, pch = 17, col = cols[1:n.series])
          text(x, y, 1:n.series, pos = 3)
        }
        plot.new()
        color.legend(0, 0, 1, 1, rect.col = col.palette(20),
            legend = format(pretty(cube.info[pars.plot[i], , , fourth.dim.t], digits = 5), scientific = -1)
            , gradient = 'x', align = 'rb', cex = 0.7)
      } else {
        plot.new()
        plot.new()
      }   
    } 
    mtext(outer = TRUE, line = 0, side = 3, text = sub('.data', '', pars.plot)
        , at = 1 / length(pars.plot) * (0:(length(pars.plot) - 1)) + 1 / (2 * length(pars.plot)))
    mtext(outer = TRUE, side = 3, line = 1, cex = 3, text = paste('dim ', fourth.dim.t))
    
    ## plot example series
    ind.series <-  array(FALSE, dim = dim(data.cube.sort)[1:2])
    ind.series[as.matrix(data.frame(a = grids.valid[ind.rand[1:n.series],1], b = grids.valid[ind.rand[1:n.series],2]))]=TRUE
    plotBG(rgb(0.9,0.9,0.9))
    data.cube.t <- data.cube.sort[,,,fourth.dim.t]
    if (cube.info.agg[fourth.dim.t, 'series.empty'] != 1) {
      ydata = t(scale(array(data.cube.t[indexDatacube(datacube = data.cube.t, logical.ind = ind.series, dims = c(1,2))],
                            dim = c(length(time), n.series)), scale =  FALSE, center =  TRUE))
      plotNLines(x.data = 1:length(time), y.data = ydata, option = 'stacked', colors = cols[1:n.series], scale = 0.8, ...)
      legend('topright', legend = paste(1:n.series, apply(ind.orig[1:n.series, ], MARGIN = 1, paste, collapse = ','), sep = ': '),
          col = cols[1:n.series], lty = lty[1:n.series], lwd = 2, bg = 'white', cex = par()$cex*2, ...)
    }
  }
  if (!(class(data.object)=='NetCDF'))
    close.nc(file.con)
  cat('Finished!\n')
  ##value<<
  ## some overview statistics of the different datacubes.
  invisible(cube.info.agg)
}

