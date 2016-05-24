decomposeNcdf = structure(function(
    ##title<< Spectrally decompose all time series in a netCDF datacube
  file.name             ##<< character: name of the ncdf file to decompose. The file has to be in the current working directory!
  , borders.wl          ##<< list: borders of the different periodicity bands to extract. Units are
                        ##   sampling frequency of the series. In case of monthly data border.wl<- list(c(11, 13))
                        ##   would extract the annual cycle (period = 12). For details, see the documentation of filterTSeriesSSA.
  , calc.parallel = TRUE##<< logical: whether to use parallel computing. Needs package doMC process.
  , center.series = TRUE##<< SSA calculation parameter: see the documentation of filterTSeriesSSA!
  , check.files = TRUE  ##<< logical: whether to use checkNcdfFile to check ncdf files for consistency.
  , debugging = FALSE   ##<< logical: if set to TRUE, debugging workspaces or dumpframes are saved at several stages
                        ##   in case of an error.     
  , harmonics = c()     ##<< SSA calculation parameter: Number of harmonics to be associated with each band. See the
                        ##   documentation of filterTSeriesSSA!
  , M = c()             ##<< SSA calculation parameter. Window length for time series embedding (can be different
                        ##   for each element in borders.wl): see the documentation of filterTSeriesSSA.
  , max.cores = 16      ##<< integer: maximum number of cores to use.
  , n.comp = c()        ##<< SSA calculation parameter: see the documentation of filterTSeriesSSA!
  , pad.series = c(0,0) ##<< SSA calculation parameter: see the documentation of filterTSeriesSSA!
  , print.status = TRUE ##<< logical: whether to print status information during the process
  , ratio.const = 0.05  ##<< numeric: max ratio of the time series that is allowed to be above tresh.const for the time series
                        ##   still to be not considered constant. 
  , repeat.extr = rep(1,times=length(borders.wl))##<< SSA calculation parameter: see the documentation of filterTSeriesSSA!
  , tresh.const = 1e-12 ##<< numeric: value below which abs(values) are assumed to be constant and excluded
                        ##   from the decomposition
  , var.names = 'auto'  ##<< character string: name of the variable to fill. If set to 'auto' (default), the name
                        ##   is taken from the file as the variable with a different name than the dimensions. An
                        ##   error is produced here in cases where more than one such variables exist.
  , ...                 ##<< additional arguments transferred to filterTSeriesSSA.
  )
  ##description<<
  ## Wrapper function to automatically decompose gridded time series inside a ncdf file and save the results
  ## to another ncdf file using SSA.
  ##
  ##details<<
  ## This is a wrapper function to automatically load, decompose and save a ncdf file using Singular Spectrum Analysis
  ## (SSA). It facilitates parallel computing and uses the filterTSeriesSSA() function. Refer to
  ## the documentation of filterTSeriesSSA() for details of the calculations and the necessary parameters, especially
  ## for how to perform stepwise filtering.
  ##
  ##
  ## NCDF file specifications
  ##
  ## Due to (possible) limitations in file size the ncdf file can only contain one variable and the one dimensional
  ## coordinate variables. The file has to contain one time dimension called 'time'. This function will
  ## create a second ncdf file identical to the input file but with an additional dimension called 'spectral.bands'
  ## which contains the separated spectral bands. In general the data is internally split into
  ## individual time series along ALL dimensions other than time, e.g. a spatiotemporal data cube would be separated
  ## into individual time series along its longitude/latitude dimension . The individual series are  decomposed
  ## and finally combined, transposed and saved in the new file.
  ##
  ## The NCDF file may contain NaN values at
  ## grid locations where no data is available (e.g. ocean tiles) but individual time series from single "valid"
  ## grid points must not contain missing values.
  ## In other words, decomposition is only performed for series without missing values, results for non gap-free series
  ## will be missing_value the results file.
  ##
  ## The function has only been exhaustively tested with ncdf files with two spatial dimensions (e.g. latitude and
  ## longitude) and the time dimension. Even though it was programmed to be more flexible, its functionality can not
  ## be guaranteed under circumstances with more and/or different dimensions.
  ## Input NCDF files should be compatible with the Climate Forcasting (CF) 1.5 ncdf conventions. Several crucial
  ## attributes and dimension units are checked and an error is caused if the convention regarding these aspects is
  ## not followed. Examples are the attributes scale_factor, add_offset _FillValue and the units for the time dimension
  ##
  ## Parallel computing
  ##
  ## If calc.parallel == TRUE, single time series are decomposed with parallel computing. This requires
  ## the package doMC  (and its dependencies) to be installed on the computer.
  ## Parallelization with other packages is theoretically possible but not yet implemented. If
  ## multiple cores are not available, setting calc.parallel to FALSE will cause the process to be
  ## calculated sequential without these dependencies. The package foreach is needed in all cases.

  ##seealso<<
  ##\code{\link[Rssa]{ssa}}, \code{\link[spectral.methods]{filterTSeriesSSA}}, \code{\link{gapfillNcdf}}

  ##value<<
  ##Nothing is returned but a ncdf file with the results is written in the working directory.

  ############################################################################################################
{
  ##TODO add mechanism to get constant values in datacube after calculation.
  ##TODO Try zero line crossings for frequency determination
  ##TODO Make method reproducible (seed etc)
  ##TODO Add way to handle non convergence
  
  ## prepare parallel back end
  if (calc.parallel) 
    w <- registerParallel('doMC', max.cores)
  
  ##save argument values of call
  args.call.filecheck <- as.list(environment())
  args.call.global    <- convertArgs2String()
  if (print.status & !interactive()) {
    print('Arguments supplied to function call:')
    print(paste(paste(names(args.call.filecheck), args.call.filecheck, sep = ':'), collapse = '; '))
  }
  
  ## check input
  res.check     <- do.call(.checkInputNcdfSSA,
                           c(SSAprocess = 'Decompose', args.call.filecheck))
  file.con.orig <- res.check$file.con.orig    
  if (length(var.names) == 1 && var.names == 'auto') 
    var.names   <- readNcdfVarName(file.name)
  
  ## open ncdf files
  if (print.status)
    cat(paste(Sys.time(), ' : Creating ncdf file for results. \n', sep=''))
  file.name.copy  <- paste(sub('.nc$','',file.name), '_specdecomp.nc', sep='')
  file.con.copy   <- create.nc(file.name.copy)
  Sys.chmod(file.name.copy, mode = "0777")

  ## set default parameters
  if (!calc.parallel)
    max.cores                   <- 1
  n.bands                       <- length(unlist(borders.wl))-length(borders.wl)
  n.steps                       <- length(borders.wl)
  n.timesteps                   <- dim.inq.nc(file.con.orig, 'time')$length
  if (length(M) == 0)
    M              <- rep(round(n.timesteps / 2,  digits = 0), times = n.steps)
  if (length(harmonics) == 0)
    harmonics      <- rep(0, times = n.steps)
  if (length(n.comp) == 0)
    n.comp         <- rep(50, times = n.steps)
  
  ## determine call settings for SSA
  args.call                     <- list(...)
  args.call[['borders.wl']]     <- borders.wl
  args.call[['M']]              <- M
  args.call[['harmonics']]      <- harmonics
  args.call[['n.comp']]         <- n.comp
  args.call[['pad.series']]     <- pad.series
  args.call[['print.stat']]     <- FALSE
  args.call[['plot.spectra']]   <- FALSE
  args.call[['center.series']]  <- center.series
  args.call[['repeat.extr']]    <- repeat.extr
  args.call[['debugging']]      <- debugging
  dim.values                    <- 1:n.bands
  borders.low                   <- rapply(borders.wl, function(x){x[-length(x)]})
  borders.up                    <- rapply(borders.wl, function(x){x[-1]})
  dim.name                      <- 'spectral_bands'
  
  ##prepare results file
  modifyNcdfCopyMetadata(file.con.orig, file.con.copy)    
  dim.def.nc(file.con.copy, 'spectral_bands', length(dim.values) )
  var.def.nc(file.con.copy, 'borders.low', 'NC_DOUBLE',  'spectral_bands')
  var.def.nc(file.con.copy, 'borders.up', 'NC_DOUBLE', 'spectral_bands')
  var.put.nc(file.con.copy, 'borders.low', borders.low)
  var.put.nc(file.con.copy, 'borders.up', borders.up)
  att.put.nc(file.con.copy, 'borders.up', 'long_name', 'NC_CHAR',
             'upper period border of spectral band')
  att.put.nc(file.con.copy, 'borders.low', 'long_name', 'NC_CHAR',
             'lower period border of spectral band')
  att.put.nc(file.con.copy, 'borders.up', 'unit', 'NC_CHAR',
             '[timesteps]')
  att.put.nc(file.con.copy, 'borders.low', 'unit', 'NC_CHAR',
             '[timesteps]')
  for (var.name.create in var.names) {
    var.def.nc(file.con.copy, var.name.create, var.inq.nc(file.con.orig, var.name.create)$type,
               c(var.inq.nc(file.con.orig, var.name.create)$dimids,
                 dim.inq.nc(file.con.copy, 'spectral_bands')$id))
    modifyNcdfCopyAtts(file.con.orig, var.name.create, var.name.create, file.con.copy)
  }
  sync.nc(file.con.copy)


  for (var.name in var.names) {
    data.all          <- var.get.nc(file.con.orig, var.name)
    if (print.status)
      printStatus(paste('Processing variable ', var.name, sep = ''))
    ##prepare parallel iteration parameters
    dims.ids.data       <- var.inq.nc(file.con.orig, var.names[1])$dimids + 1   
    dims.info           <- infoNcdfDims(file.con.orig)[dims.ids.data, ]
    drop.dim            <- FALSE
    if (length(dim(data.all)) > 1 ) {
      dims.cycle.id     <- sort(setdiff(1:length(dims.ids.data), match('time', dims.info$name) ))
      dims.cycle.name   <- dims.info[dims.cycle.id, 'name']
      dims.process.id   <- match('time', dims.info$name)
    } else {
      dims.cycle.id     <- 1
      dims.process.id   <- 2
      dims.cycle.name   <- 'series'
      data.all          <- array(data.all, dim = c(1, length(data.all)))
      drop.dim          <- TRUE
    }
    dims.cycle.length   <- dim(data.all)[dims.cycle.id] 
    slices.n            <- prod(dims.cycle.length)
    dims.cycle.n        <- length(dims.cycle.id)
    n.timesteps         <- dims.info[match('time', dims.info$name), 3]

    ## determine slices to process
    results.identify <- .identifyValidCellsSSA(dims.cycle.id = dims.cycle.id, dims.process.id = dims.process.id,
                                              datacube = data.all, ratio.const = ratio.const, 
                                              tresh.const = tresh.const , print.status = print.status, 
                                              slices.n = slices.n, algorithm = 'Decompose')
    iters.n = results.identify$iters.n
    slices.process = results.identify$slices.process
    values.constant = results.identify$values.constant
    slices.constant = results.identify$slices.constant 
    slices.without.gaps = results.identify$slices.without.gaps
    slices.excluded = results.identify$slices.excluded
    slices.too.gappy = results.identify$slices.too.gappy
 
    if (sum(results.identify$slices.process) == 0) {
      printStatus(paste('Specdecomp for ', var.name,' not possible. Next step ...', sep = ''))
      next
    }

    ## create 'iterator'
    args.expand.grid        <- alist()
    for (i in 1:dims.cycle.n)
      args.expand.grid[[i]] <- 1:dims.cycle.length[i]
    iter.grid.all           <- as.matrix(do.call("expand.grid", args.expand.grid))
    n.slices                <- dim(iter.grid.all)[1]
    n.iters                 <- sum(slices.process)
    max.cores               <- min(c(max.cores, n.iters))
    iter.grid               <- matrix(1, nrow = n.iters, ncol = length(dims.cycle.id) + 1)
    colnames(iter.grid)     <- c('iter.nr', dims.cycle.name)
    iter.grid[,'iter.nr']   <- 1:n.iters
    iter.grid[, dims.cycle.id + 1] <- iter.grid.all[slices.process, ]
    iters.per.cyc           <- rep(floor(n.iters/max.cores), times = max.cores)
    if (!(n.iters %% max.cores) == 0)
      iters.per.cyc[1:((n.iters %% max.cores))] <- floor(n.iters / max.cores) + 1
    iter.gridind            <- matrix(NA, ncol = 2, nrow = max.cores)
    colnames(iter.gridind)  <- c('first.iter','last.iter')
    if (max.cores > 1)  {
      iter.gridind[, 1]     <- (cumsum(iters.per.cyc) + 1) - iters.per.cyc
      iter.gridind[, 2]     <- cumsum(iters.per.cyc)
    } else {
      iter.gridind[1, ]     <- c(1, n.iters)
    }

    ## define process during iteration
    
    abind.mod=function(...)abind(..., along=3)

    ## perform calculation
    if (print.status)
      cat(paste(Sys.time(), ' : Starting calculation: Decomposing ', sum(slices.process),
                ' timeseries of length ', n.timesteps, '. \n', sep=''))
    if (calc.parallel) {
      data.results.valid.cells <- foreach(i = 1:max.cores, .combine = abind.mod, .multicombine = TRUE) %dopar% .decomposeNcdfCoreprocess(
                                                                                   iter.nr = i, print.status = print.status, data.all = data.all, 
                                                                                   n.timesteps = n.timesteps, file.name = file.name,
                                                                                   n.bands = n.bands, dims.cycle.n = dims.cycle.n,
                                                                                   iter.grid = iter.grid, args.call = args.call,
                                                                                   var.name = var.name, iter.gridind = iter.gridind, 
                                                                                   dims.process.id = dims.process.id, dims.cycle.id = dims.cycle.id)
    } else {
      data.results.valid.cells <- foreach(i = 1:max.cores
                                          ,  .combine = abind.mod,  .multicombine = TRUE) %do% .decomposeNcdfCoreprocess(iter.nr = i,
                                                                      n.timesteps = n.timesteps,  print.status = print.status,
                                                                      n.bands = n.bands, dims.cycle.n = dims.cycle.n, data.all = data.all, 
                                                                      iter.grid = iter.grid, args.call = args.call, file.name = file.name,
                                                                      var.name = var.name, dims.process.id = dims.process.id,
                                                                      dims.cycle.id = dims.cycle.id, iter.gridind = iter.gridind)
    }

    if (print.status)
      cat(paste(Sys.time(), ' : Transposing results. \n', sep=''))
    data.results.all.cells                     <- array(NA, dim = c(n.timesteps, n.bands, n.slices))
    data.results.all.cells[, , slices.process] <- data.results.valid.cells
    data.results.all.cells.trans               <- aperm(data.results.all.cells, perm = c(3, 1, 2))
    data.results.final                         <- array(as.vector(data.results.all.cells.trans),
                                                        dim = c(dims.cycle.length, n.timesteps, n.bands))
    aperm.array                                <- c(order(c(dims.cycle.id, dims.process.id)), length(c(dims.cycle.id, dims.process.id)) + 1)
    data.results.final                         <- aperm(data.results.final, aperm.array)

    if (print.status)
      cat(paste(Sys.time(), ' : Writing results to file. \n', sep=''))

    ## add missing value attribute
    if (sum(is.na(data.results.final)) > 0) {
      att.missing <- c('missing_value',  '_FillValue')[is.na(match(c('missing_value',  '_FillValue'), infoNcdfAtts(file.con.copy, var.name)[,'name']))]
      for (att.add in att.missing) {
        other.att <-  c('missing_value',  '_FillValue')[ - match(att.add, c('missing_value',  '_FillValue'))]
        if (is.element(other.att, infoNcdfAtts(file.con.copy, var.name)[,'name'])) {
          value.use <- att.get.nc(file.con.copy,  var.name,  other.att)
        } else {
          value.use <- -9999.0
        }
        att.put.nc(file.con.copy, var.name, att.add, var.inq.nc(file.con.copy, var.name)$type, value.use )
      }
    }

    ## save results
    if (drop.dim)
      data.results.final <- drop(data.results.final)
    var.put.nc(file.con.copy, var.name, data.results.final)
    sync.nc(file.con.copy)
  }
  
  ## add attributes with process information to ncdf files
  all.args     <- formals(filterTSeriesSSA)
  all.args[na.omit(match(names(args.call), names(all.args)))] <- args.call[is.element(names(args.call), names(all.args))]
  red.args     <- all.args[c('borders.wl', 'M', 'n.comp', 'harmonics', 'tolerance.harmonics',
                             'var.thresh.ratio', 'grouping', 'pad.series', 'SSA.methods', 'repeat.extr')]
  string.args  <- paste(paste(names(red.args), sapply(red.args, function(x)paste(x, collapse=', '))
                              , sep=': '), collapse='; ')
  att.put.nc(file.con.copy, 'NC_GLOBAL', 'Decomposed_by', 'NC_CHAR', getSysinfo())
  att.put.nc(file.con.copy, 'NC_GLOBAL', 'Decomposed_on', 'NC_CHAR', as.character(Sys.time()))
  att.put.nc(file.con.copy, 'NC_GLOBAL', 'Decomposition_settings', 'NC_CHAR', string.args)
  hist.string.append <- paste('Spectrally decomposed on ', as.character(Sys.time()),
                              ' by ', Sys.info()['user'], sep='')
  if (is.element('history', infoNcdfAtts(file.con.copy, 'NC_GLOBAL')[, 'name'])) {
    hist.string <- paste(att.get.nc(file.con.copy, 'NC_GLOBAL', 'history'), '; ', hist.string.append)
    att.put.nc(file.con.copy, 'NC_GLOBAL', 'history', 'NC_CHAR', hist.string)
  } else {
    att.put.nc(file.con.copy, 'NC_GLOBAL', 'history', 'NC_CHAR', hist.string.append)
  }
  close.nc(file.con.copy)
  close.nc(file.con.orig)
  if (print.status)
    cat(paste(Sys.time(), ' : Calculation successfully finished. \n', sep=''))
  return(list(finished =  TRUE))
}
, ex = function(){
  ## Example for the filtering of monthly data
  filename   <- '<filename>.nc'
  # Extract yearly cycle, intra annual part and high frequency residual in several steps
  borders.wl <- list(a = c(10, 14)
                     , b = c(12, Inf)
                     , c = c(0, 12))
  M         <- c(2*12, 4*12, 12)
  #extract first four harmonics for yearly cycle
  harmonics <- c(4, 0, 0)

  # uncomment and run 
  # decomposeNcdf(file.name = filename, borders.wl = borders.wl, M = M, harmonics = harmonics)

  # Extract yearly cycle, intra annual part and high frequency residual in one step
  borders.wl <- list(c(0,10,14,Inf))
  # use the same M for all bands
  M          <- c(2*12)
  # uncomment and run
  #decomposeNcdf(file.name = filename, borders.wl = borders.wl, M = M)
})

#####################################      core process       #############################

.decomposeNcdfCoreprocess = function(iter.nr, n.timesteps, n.bands, dims.cycle.n, data.all, 
                                    iter.grid, args.call, var.name, dims.process.id,
                                    dims.cycle.id, iter.gridind, print.status, file.name)
  ##description<<
  ## Helper function for decomposeNcdf. Is run parallellised on each core.
{
  iter.ind                   <- iter.gridind[iter.nr, ]
  data.results.iter          <- array(NA,dim=c(n.timesteps, n.bands, diff(iter.ind) + 1))
  for (j in 1:(diff(iter.ind) + 1)) {
    ind.total = (iter.ind[1]:iter.ind[2])[j]
    data.results.iter.t=try(
      {
        ind.extract <- list(data.all)
        for (i in 1:length(dims.cycle.id))
          ind.extract[[dims.cycle.id[i] + 1]] <- iter.grid[ind.total, i + 1]
        ind.extract[[dims.process.id + 1]] <- TRUE
        args.call.t             <- args.call
        args.call.t[['series']] <- as.numeric(do.call('[', ind.extract))
        series.decomp           <- do.call(filterTSeriesSSA, args.call.t)
        t(series.decomp$dec.series)
      })
    if (class(data.results.iter.t) == 'try-error') {
        print(paste('Error occoured at iteration ', iter.nr, ' and loop ', j, '!', sep=''))
        error.from.calc                 <- data.results.iter.t
        trace.save                      <- traceback()
        error.from.calc                 <- data.results.iter.t
        data.results.iter.t             <- matrix(Inf, ncol=n.bands, nrow=n.timesteps)
        system.info                     <- sessionInfo()
        path.file                       <- file.path('/', 'Net', 'Groups', 'BGI', 'tmp', 
                                                     'jbuttlar', 'Cluster_jobs_debugging', sub('/Net/Groups/BGI/', '', getwd()))
        if (!file.exists(path.file))
          system(paste('mkdir -p ', path.file, sep = ''))                      
        file.name.t                     <- file.path(path.file, paste('workspace_error_', file.name, '_',
                                                                      iter.nr, '_', j, sep = ''))
        print(paste('Saving workspace to file ', file.name.t, '.rda', sep = ''))
        dump.frames(to.file = TRUE, dumpto = file.name.t)
      }
    data.results.iter[,,j]  <- data.results.iter.t
    if (print.status)
      if (iter.nr == 1 &&( diff(iter.ind) < 20  || (j%%(ceiling((diff(iter.ind)) / 20)) == 0)))
        cat(paste(Sys.time(), ' : Finished ~', round(j / (diff(iter.ind) + 1) * 100), '%. \n', sep=''))
  }
  return(data.results.iter)
}
globalVariables('i')
