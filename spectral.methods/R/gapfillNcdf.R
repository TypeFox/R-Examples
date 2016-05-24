gapfillNcdf <- structure(function(
##title<< Fill gaps in time series or spatial fields inside a netCDF file using SSA.
##description<< Wrapper function to automatically fill gaps in series or spatial fields
##              inside a ncdf file and save the results to another ncdf file. This can
##              be done via 1D,  2D or the spatio - temporal gap 3D filling algorithm
##              of Buttlar et. al (2014).
amnt.artgaps = rep(list(   rep(list(c(0.05, 0.05)), times = length(dimensions[[1]]))) , times = length(dimensions)) ##<< list of numeric vectors: 
                      ##             The relative ratio (length gaps/series length) of
                      ##             artificial gaps to include in the "innermost" SSA loop (e.g. the value used by the
                      ##             SSA run for each individual series/slice). These ratio is used to determine 
                      ##             the iteration with the best prediction (c(ratio big gaps, ratio small gaps)) (see ?gapfillSSA for details )                                  
, amnt.iters = rep(list(   rep(list(c(10, 10)), times = length(dimensions[[1]]))) , times = length(dimensions)) ##<< list of integer vectors: 
                      ##             amount of iterations performed for the outer and inner
                      ##             loop (c(outer,inner)) (see ?gapfillSSA for details)
, amnt.iters.start = rep(list(   rep(list(c(1, 1)), times = length(dimensions[[1]]))) , times = length(dimensions)) ##<< list of integer vectors:
                      ##             index of the iteration to start with (outer, inner). If this
                      ##             value is >1, the reconstruction (!) part is started with this iteration. Currently
                      ##             it is only possible to set this to values >1 if amnt.artgaps and size.biggap == 0.
, calc.parallel = TRUE##<< logical:  whether to use parallel computing. Needs the packages doMC and foreach
                      ##             to be installed.                                  
, debugging = FALSE   ##<< logical: if set to TRUE, debugging workspaces or dumpframes are saved at several stages
                      ##            in case of an error.
,  debugging.SSA =  FALSE
, dimensions = list(list('time'))    ##<< list of string vectors: 
                      ##             setting along which dimensions to perform SSA. These names
                      ##             have to be identical to the names of the dimensions in the ncdf file. Set this to
                      ##             'time' to do only temporal gap filling or to (for example) c('latitude','longitude')
                      ##             to do 2 dimensional spatial gap filling. See the description for details on how to
                      ##             perform spatio-temporal gap filling.
, file.name           ##<< character: name of the ncdf file to decompose.  The file has to be in the current working directory!                     
, first.guess = 'mean'##<< character string: if 'mean', standard SSA procedure is followed (using zero as the first guess). 
                      ##             Otherwise this is the file name of a ncdf file with the same dimensions
                      ##             (with identical size!) as the file to fill which contains values used as a
                      ##             first guess (for the first step of the process!). This name needs to be exactly
                      ##             "<filename>_first_guess_step_1.nc".                                  
, force.all.dims = FALSE ##<< logical: if this is set to true, results from dimensions not chosen as the best guess are used
                      ##             to fill gaps that could not be filled by the best guess dimensions due to too gappy slices etc. .                                  
, gaps.cv = 0         ##<< numeric: ratio (between 0 and 1) of artificial gaps to be included prior to the first
                      ##             cross validation loop that are used for cross validation.
, keep.steps = TRUE   ##<< logical: whether to keep the files with the results from the single steps                                                        
, M                   ##<< list of single integers: window length  or embedding dimension in time steps. If not
                      ##             given,  a default value of 0.5*length(time series) is computed.                                                               
, max.cores = 8       ##<< integer: maximum number of cores to use (if calc.parallel = TRUE).                                  
, max.steps = 10      ##<< integer: maximum amount of steps in the variances scheme                                 
, n.comp = lapply(amnt.iters, FUN = function(x)x[[1]][[1]][1] * 2) ##<< list of single integers: 
                      ##             amount of eigentriples to extract (default (if no values are
                      ##             supplied) is 2*amnt.iters[1] (see ?gapfillSSA for details)
, ocean.mask = c()    ##<< logical matrix: contains TRUE for every ocean grid cell and FALSE for land cells. Ocean grid
                      ##             cells will be set to missing after spatial SSA and will be excluded from temporal SSA.
                      ##             The matrix needs to have dimensions identical in order and size to the spatial dimensions
                      ##             in the ncdf file. As an alternative to a R matrix, the name of a ncdf file can be supplied.
                      ##             It should only contain one non coordinate variable with 1 at ocean cells and 0 at land cells.
, pad.series = rep(list(   rep(list(c(0, 0)), times = length(dimensions[[1]]))) , times = length(dimensions)) ##<< list of integer vectors (of size 2):
                      ##             length of the extracts from series to use for
                      ##             padding. Only possible in the one dimensional case. See the documentation of gapfillSSA for details!
, print.status = TRUE ##<< logical: whether to print status information during the process                                  
, process.cells = c('gappy','all')[1] ##<< character string: 
                      ##             which grid/series to process. 'gappy' means that only series grids with actual
                      ##             gaps will be processed, 'all' would result in also non gappy grids to be subjected to SSA. The
                      ##             first option results in faster computation times as reconstructions for non gappy grids/series
                      ##             are technically not needed for gap filling, whereas the second option provides a better
                      ##             understanding of the results trajectory to the final results.                                  
, process.type = c('stepwise', 'variances')[1] ##character string<< Use 'stepwise' to stepwise fill using one dimension setting
                      ##             and 'variances' to alternate between different dimension settings using the one with
                      ##             the lowest residual variance to proceed (e.g. for the so called spatiotemporal 3D scheme). 
, ratio.const = 0.05  ##<< numeric: max ratio of the time series that is allowed to be above tresh.const for the time series
                      ##             still to be not considered constant.                                 
, ratio.test = 1      ##<< numeric: ratio (0-1) of the data that should be used in the cross validation step. If set to 1,
                      ##             all data is used.
, reproducible = FALSE##<< logical: Whether a seed based on the characters of the file name should be set
                      ##            which forces all random steps, including the nutrlan SSA algorithm to be
                      ##            exactly reproducible.
, size.biggap = rep(list(   rep(list(20) , times = length(dimensions[[1]]))) , times = length(dimensions)) ##<< list of single integers: 
                      ##             length of the big artificial gaps (in time steps) (see ?gapfillSSA for details)
, tresh.const = 1e-12 ##<< numeric: value below which abs(values) are assumed to be constant and excluded
                      ##             from the decomposition.
, tresh.converged = 0 ##<< numeric: ratio (0-1): determines the amount of SSA iterations that have to converge so that no error
                      ##             is produced.
, tresh.fill = 0.1     ##<< numeric fraction (0-1):
                      ##             This value determines the fraction of valid values below which
                      ##             series/grids will not be filled in this step and are filled with the first guess from the
                      ##             previous step (if any). For filling global maps while using a ocean.mask you need
                      ##             to set this value relative to the global grid size (and not only the land mask). Setting this
                      ##             value to zero would mean that also slices/series without any "real" values are tried to
                      ##             be filled with the "first guess" value of the last iteration alone. This can only be done
                      ##             if the cross validation scheme is switched off (e.g. by setting amnt.artgaps and size.biggap
                      ##             to zero.
, var.names = 'auto'  ##<< character string: name of the variable to fill. If set to 'auto' (default), the name
                      ##             is taken from the file as the variable with a different name than the dimensions. An
                      ##             error is produced here in cases where more than one such variable exists.
, ...                 ##<< further settings to be passed to gapfillSSA
)
##details<<
## This is a wrapper function to automatically load, gapfill and save a ncdf file using SSA.
## Basically it automatically runs gapfillSSA (see also its documentation) for all                         
## time series or grids in a ncdf file automatically. Theoretically and  
## conceptionally all methods could also be applied to simple datacubes (i.e. R arrays) 
## and not only ncdf files. However,  this has not yet been implemented. The values
## for several function arguments have to be supplied as rather complicated nested lists                         
## to facilitate the usage of different values per step (see 'stepwise calculation' for details)
##                         
## dimensions:
## It is generally possible to perform one,  two,  or spatiotemporal 3D SSA
## (as in Buttlar et al (2014)) for gap filling.  This is set by using the argument
## 'dimensions'. If only one string corresponding to a dimension name in the target
## ncdf file is supplied, only  vectors in the direction of this dimension are
## extracted and filled. If two dimension names are supplied,
## matrices  (i.e. spatial grids) along these dimension are extracted and 2D SSA
## is computed. To start the 3D spatio - temporal scheme (Buttlar et.al (2014)) which
## iterates between 1D temporal and 2D spatial SSA,  set
## dimensions = list(list("time", c("longitude","latitude"))). 
##
## stepwise calculation:
## The algorithm can be run step wise with different settings for each step
## where the results from each step can be used as 'first guesses' for the subsequent
## step. To do this, amnt.artgaps, size.biggap, amnt.iters, n.comp, M, tresh.fill
## and dimensions have to be supplied as lists. For each nth iteration step the
## values of the corresponding nth list element will be used. At each of these
## nth iteration steps, several repetitions with different dimensions are possible
## (as is the case with the 3D spatiotemporal scheme). To facilitate this, the 
## individual list elements at each step have to be lists containing the different
## dimension names. As a consequence, all these arguments have to be nested lists.
## This is the case also if only one dimension is used (i.e. to do only temporal
## 1D SSA,  dimensions = list(list('time')).
## One example for an application where supplying different settings for all these
## steps would be user defined spatio-temporal gap filling. This allows to clearly
## define which dimension (and M, gap amount etc) to use at each step of the
## process. On the contrary,  the spatiotemporal gapfilling method applied by
## Buttlar et. al. (2014) uses identical settings for each (outer loop)
## iteration step and automatically determines which dimension to use. For
## this procedure the first list element of dimensions (and the other stepwise
## arguments) is recycled during each step. Hence, a list of only length one has
## to be supplied (dimensions = list(list(c('longitude','latitude'),'time')) )
## (see details for dimensions above).
##
##
## NCDF file specifications:
## Due to limitations in the file size, the ncdf file has to contain one variable (and the dimensional
## coordinate variables) (for the time being). This function will
## create a second ncdf file identical to the input file but with an additional variable called 'flag.orig',
## which contains zero for gapfilled values and 1 for not filled values.
## The function has only been tested with ncdf files with two spatial dimensions (e.g. lat and long) and
## one time dimension. Even though it was programmed to be more flexible, its functionality can not
## be guaranteed under circumstances with more dimensions.                         
##
## non fillable gid cells:
## Using the 3D method would result in a completely filled datacube. To prevent
## the filling of grid cells where no reasonable guess via the gapfilling may be
## achieved (i.e. ocean grid cells in the case of terrestrial data), a matrix
## indicating these grid cells can be supplied (see 'ocean.mask')
##

##
## Parallel computing:
## If calc.parallel==TRUE, single time series are filled with parallel computing. This requires
## the package doMC (and its dependencies) to be installed on the computer.
## Parallelization with other packages is theoretically possible but not yet implemented. If
## multiple cores are not available, setting calc.parallel to FALSE will cause the process to be
## calculated sequential without these dependencies. The package foreach is needed in all cases.


##references<<
## v. Buttlar, J., Zscheischler, J., and Mahecha, M. D. (2014): An extended approach for
## spatiotemporal gapfilling: dealing with large and systematic gaps in geoscientific
## datasets, Nonlin. Processes Geophys., 21, 203-215, doi:10.5194/npg-21-203-2014                       
##seealso<<
##\code{\link[Rssa]{ssa}}, \code{\link[spectral.methods]{gapfillSSA}}, \code{\link{decomposeNcdf}}

                 
##value<<
##Nothing is returned but a ncdf file with the results is written.
######################################################################################################
{
    ##TODO remove aperm steps
    ##TODO extract iloop convergence information for all loops
    ##TODO test inner loop convergence scheme for scenarios 
    ##TODO indicate fraction of gaps for each time series
    ##TODO break down world into blocks
    ##TODO integrate onlytime into one dimensional variances scheme
    ##TODO facilitate one step filling process with global RMSE calculation
    ##TODO save convergence information in ncdf files
    ##TODO check for too gappy series at single dimension setting
    ##TODO create possibility for non convergence and indicate this in results
    ##TODO facilitate run without cross validation repetition
    ##TODO test stuff with different dimension orders in the file and in settings
    ##TODO substitute all length(processes)==2 tests with something more intuitive
    ##TODO put understandable documentation to if clauses
    ##TODO remove first guess stuff
    ##TODO incorporate non convergence information in final datacube
    ##TODO facilitate easy run of different settings (e.g. with different default settings)
    ##TODO switch off "force.all.dims" in case of non necessity
    ##TODO delete/modify MSSA stuff

    ##obsolete MSSA stuff  
    MSSA =  rep(list(   rep(list(FALSE), times = length(dimensions[[1]]))) , times = length(dimensions)) 
    MSSA.blck.trsh = tresh.fill                                    
    MSSA.blocksize = 1

    #save argument values of call
    printStatus(paste('Filling file', file.name))
    args.call.filecheck <- as.list(environment())
    args.call.global    <- convertArgs2String()
    if (print.status & !interactive()) {
      print('Arguments supplied to function call:')
      cat(paste(paste(paste(names(args.call.filecheck), args.call.filecheck, sep=':'),
                      collapse = '; '), '\n\n'))
    }
    
    #set seed based on file name
    if (reproducible) {
       file.name.cl <-  gsub('[[:punct:]]', '', file.name)
       ind.rev      <-  round(seq(1, nchar(file.name.cl), length.out = 5), digits = 0)
       seed.big     <-  as.numeric(paste(match(unlist(strsplit(file.name.cl, ''))[ind.rev], 
                                   c(letters, LETTERS, 0:9)) , collapse = '' )   )
       seed         <-  (seed.big)%%(2^32)
       if (print.status)
         cat(paste(Sys.time(), ' : Using seed ', seed, ' for calculations.\n', sep=''))
       set.seed(seed)
    } else {
      seed <- c()
    }

    # necessary variables
    if (gaps.cv != 0) {
      processes <- c('cv', 'final')
    } else if (gaps.cv == 0) {
      processes <- c('final') 
      art.gaps.values            <- NULL            
    }
    process.converged.detail = list()
    if (process.type == 'variances') {
      step.chosen         <- matrix(NA, 2, max.steps)
      iter.chosen         <- array(NA, dim=c(3, max.steps, length(processes)))
      process.converged   <- matrix(NA, length(processes), max.steps)      
      n.steps             <- max.steps
      var.steps           <- array(NA, dim = c(max.steps, max(unlist(n.comp)), 
                                         length(dimensions[[1]])))
    } else {
      n.steps             <- length(dimensions)
      step.chosen         <- matrix(NA, 2, n.steps)
      iter.chosen         <- array(NA, dim=c(3, n.steps, length(processes)))
      process.converged   <- matrix(NA, length(processes), n.steps)      
      step.chosen[1, ]    <- 1
      step.chosen[2, ]    <- 1:n.steps
    }
    dimnames(step.chosen) <- list(c('dim','step'), paste('step', 1:dim(step.chosen)[2]))
    dimnames(iter.chosen) <- list(c('outer','inner', 'amnt.na'),
                                  paste('step', 1:dim(iter.chosen)[2]), processes)
    dimnames(process.converged) <- list(processes)


    # debugging and information variables
    finished              <- FALSE
    iterpath              <- data.frame(time = Sys.time(), var.name = 'none',
                                        process = 'none', step = 0, calc.repeat = 0,
                                        dimensions = 0,  otherdim = NA, fg.filename.nxt = '',
                                        fg.used = FALSE, stringsAsFactors = FALSE)
    included.otherdim     <- rep(FALSE, n.steps)
    args2SSA              <- list()

    
    #check input, check first guess, transfer and check ocean mask
    res.check     <- do.call(.checkInputNcdfSSA,
                             c(SSAprocess = 'Gapfill', args.call.filecheck))
    ocean.mask    <- res.check$ocean.mask

    #determine some filenames etc
    if (var.names[1] == 'auto')
      var.names = readNcdfVarName(file.name)
    file.con.orig <- open.nc(file.name)
    dims.info     <- infoNcdfDims(file.con.orig)[pmatch(c('lat', 'lon', 'time'), infoNcdfDims(file.con.orig)$name), ]
    dims.info$id <- 0:(dim(dims.info)[1] - 1)

    
    ## start parallel processing workers
    if (calc.parallel) 
      w <- registerParallel('doMC', max.cores)

    # start loops
    for (var.name in var.names) {
      result.rotation <-  transNcdfRotate(data.object = file.con.orig,
                                          var.name =  var.name,  reverse.dim =  FALSE)
      datacube        <-  result.rotation$data
      aperm.reverse   <-  result.rotation$aperm.reverse
      if (sum(is.na(datacube)) == 0) {
        print(paste(var.name, 'does not contain gaps. Gap filling for this data omitted!'))
        next
      } 
      ind.artgaps.out  <- array(FALSE, dim = dim(datacube))
      printStatus(paste('Filling variable ', var.name, sep = ''))
      for (process in processes) {
        
        ## insert gaps for cross validation      
        if (process == 'cv') {
          if (print.status)
            cat(paste(Sys.time(), ' : Starting cross validation loop. \n', sep = ''))        
          indices.t                  <- sample(which(!is.na(datacube)),
                                               floor(gaps.cv * sum(!is.na(datacube))))
          ind.artgaps.out[indices.t] <- TRUE
          art.gaps.values            <- datacube[ind.artgaps.out]
          datacube[ind.artgaps.out]  <- NA
        } else if (process == 'final') {
          if (print.status)
            cat(paste(Sys.time(), ' : Starting final filling loop. \n', sep = ''))
          datacube[ind.artgaps.out]  <- art.gaps.values
          if (length(processes) == 2)
            n.steps  <- max(step.chosen['step',])
        }
        for (h in 1:n.steps) {
          if (print.status)
            cat(paste(Sys.time(), ' : Starting step ', h, '\n',sep = ''))
          
          ## determine different iteration control parameters
          if (process.type == 'stepwise') {
            ind                   <- h     
            n.dims.loop           <- length(dimensions[[ind]])
            n.calc.repeat         <- 1
            dims.calc             <- 1:n.dims.loop
            ratio.test.t          <- 1
          } else if (process.type == 'variances') {
            ind                   <- 1
            n.dims.loop           <- length(dimensions[[ind]])
            if (process == 'final') {
              if (force.all.dims) {
                dims.calc         <- 1:n.dims.loop
              } else if (!force.all.dims) {  
                if (length(processes) == 2) {
                  dims.calc       <- step.chosen['dim', h]
                } else if (length(processes) == 1) {
                  dims.calc       <- 1
                }
              }
              n.calc.repeat       <- 1
              ratio.test.t        <- 1
            } else if (process == 'cv') {
              dims.calc           <- 1:n.dims.loop
              if (force.all.dims) {
                ratio.test.t      <- 1
              } else {
                ratio.test.t      <- ratio.test 
              }
            }                  
            if (ratio.test.t == 1) {
              n.calc.repeat       <- 1
            } else {
              n.calc.repeat       <- 2
            }	        
          }
          if (!exists('pred.measures')) {
            pred.measures <- array(NA, dim = c(3, n.steps, n.dims.loop))
            dimnames(pred.measures) <- list(c('var.res.steps','RMSE','MEF'), 
                                            paste('step',1:n.steps), 
                                            sapply(args.call.filecheck$dimensions[[1]],
                                                   function(x) paste(x, collapse='-')))
          }
          
          for (g in 1:n.calc.repeat) {
            for (l in dims.calc) {
              if (g == 2 && step.chosen['dim', h] != l)
                next
              tresh.fill.dc         <- tresh.fill
              
              ##prepare parallel iteration parameters
              if (process.type == 'stepwise') {
                amnt.iters.loop       <- amnt.iters[[h]][[1]]
                amnt.iters.start.loop <- amnt.iters.start[[h]][[1]]
              } else if (process.type == 'variances') {
                if (h == 1)
                  tresh.fill.dc       <- tresh.fill
                amnt.iters.loop       <- c(h, amnt.iters[[1]][[l]][2])
                amnt.iters.start.loop <- c(h, 1)
              }
              if (print.status)
                cat(paste(Sys.time(), ' : Starting process for filling dimension: ',
                          paste(dimensions[[ind]][[l]], collapse=','), ' \n', sep = ''))
              iterpath <- rbind(iterpath, data.frame(time = Sys.time(), var.name = var.name
                                                     , process = process, step = h, calc.repeat = g
                                                     , dimensions = paste(dimensions[[ind]][[l]],
                                                         collapse = ','),
                                                     otherdim = NA,  fg.filename.nxt= NA, fg.used = FALSE))
              drop.dim = FALSE
              dims.process        <- dimensions[[ind]][[l]]
              dims.process.id     <- dims.info[match(dims.process, dims.info$name), 1]
              dims.process.length <- dims.info[match(dims.process, dims.info$name), 3]
              if (length(dim(datacube)) == 1 &&  dims.process == 'time') {
                dims.cycle.id     <- 0
                dims.cycle.amnt   <- 1
                dims.cycle.length <- 1
                dims.process.id   <- 1
                dims.cycle        <- 'series'
                datacube          <- array(datacube, dim = c(1, length(datacube)))
                drop.dim          <- TRUE
              } else {
                dims.cycle        <- dims.info[ - match(dims.process, dims.info$name), 2]
                dims.cycle.id     <- dims.info[match(dims.cycle, dims.info$name), 1]
                dims.cycle.amnt   <- length(dims.cycle.id)
                dims.cycle.length <- dims.info[match(dims.cycle, dims.info$name), 3]
              }  
              datapts.n           <- prod(dims.info[match(dims.process, dims.info$name), 3])
              slices.n            <- prod(dims.cycle.length)

              
              diminfo.step <- list(dims.process = dims.process,
                                   dims.process.id = dims.process.id,
                                   dims.process.length = dims.process.length,
                                   dims.cycle.id = dims.cycle.id,
                                   dims.cycle.length = dims.cycle.length)
              ##determine call settings for SSA
              iterinf = list(process = process, h = h, g = g, l = l)
              args.call.SSA <- list(amnt.artgaps = amnt.artgaps[[ind]][[l]],
                                    M = M[[ind]][[l]],
                                    size.biggap = size.biggap[[ind]][[l]],
                                    n.comp = n.comp[[ind]][[l]],
                                    pad.series = pad.series[[ind]][[l]],
                                    amnt.iters = amnt.iters.loop,
                                    amnt.iters.start = amnt.iters.start.loop,
                                    print.stat   = FALSE,
                                    plot.results = FALSE,
                                    debugging = debugging.SSA, seed = seed,
                                    iterinf = iterinf)
              data.step <- list(iterinf = iterinf, args = args.call.SSA)
              args2SSA[[length(args2SSA) + 1]] <- data.step
              
              
              ##get first guess
              
              if (h > 1 && exists('file.name.guess.next')) {
                file.con.guess   <- open.nc(file.name.guess.next)
                first.guess      <- transNcdfRotate(data.object = file.con.guess,
                                                    var.name = var.name, reverse.dim =  FALSE)$data
                close.nc(file.con.guess)
                iterpath[dim(iterpath)[1], 'fg.filename.nxt'] <- file.name.guess.next
                iterpath[dim(iterpath)[1], 'fg.used'] <- TRUE            
              }
              ##run calculation
              ##TODO try to stop foreach loop at first error message!
              args.Datacube <- c(list(datacube = datacube, max.cores = max.cores,
                                      tresh.fill.dc = tresh.fill.dc, first.guess = first.guess,
                                      dims.process.id = dims.process.id, dims.cycle.id = dims.cycle.id,
                                      dims.process = dims.process, dims.cycle = dims.cycle, 
                                      print.status = print.status, datapts.n = datapts.n, dims.info = dims.info,
                                      calc.parallel = calc.parallel, ocean.mask = ocean.mask, 
                                      debugging = debugging, h = h, l = l,  MSSA = MSSA[[ind]][[l]],
                                      MSSA.blocksize = MSSA.blocksize, ratio.test.t = ratio.test.t, g = g,
                                      MSSA.blck.trsh = MSSA.blck.trsh, file.name = file.name), 
                                 list(args.call.SSA = args.call.SSA), ratio.const = ratio.const,
                                 tresh.const = tresh.const, reproducible = reproducible)
              if (g > 1) {
                args.Datacube <- c(args.Datacube, list(slices.process = slices.process,
                                                       slices.constant = slices.constant,
                                                       values.constant = values.constant, 
                                                       slices.excluded = slices.excluded,
                                                       slices.without.gaps = slices.without.gaps))
              }
              gapfill.results   <- do.call(.gapfillNcdfDatacube, args.Datacube)
              
              gapfill.results   <- c(gapfill.results, diminfo.step)
              if (is.null(gapfill.results$reconstruction) && is.null(gapfill.results$data.variances) &&
                  n.steps == 1) {
                print(paste('No series available for filling and only one step process',
                            'chosen. Stopping gapfilling.',  sep = ''))
                return(list(finished = FALSE))
              }  
              if (process.type == 'variances')
                assign(paste('gapfill.results.dim', l, sep=''), gapfill.results)
            }
            

            ##test which dimension to be used for the next step
            ##TODO whole step can be excluded for "one step" processes
            if (process.type == 'variances' & ((length(processes) == 2 && process == 'cv') |
                  (length(processes) == 1 && process == 'final') )) {
              if (g == 1) {
                for (k in 1:n.dims.loop) {
                  recstr.t      <- get(paste('gapfill.results.dim', k, sep = ''))[['reconstruction']]
                  if (!is.null(recstr.t)) {
                    pred.per.t    <- var(art.gaps.values - recstr.t[ind.artgaps.out], na.rm = TRUE)
                    if (length(processes) == 2) {
                      pred.measures['var.res.steps',h ,k ] <- pred.per.t
                    } else if (length(processes) == 1) {                     # if single step process with "inner" cross validation
                      pred.measures['var.res.steps',h ,k ] <- 0
                    }  
                    if (process == 'cv') {
                      ind.test <- !is.na(recstr.t[ind.artgaps.out])
                      pred.measures['RMSE',h ,k ] <-  RMSE(art.gaps.values[ind.test],
                                                           recstr.t[ind.artgaps.out][ind.test])
                      pred.measures['MEF',h ,k ]  <-  MEF(art.gaps.values[ind.test],
                                                          recstr.t[ind.artgaps.out][ind.test])
                    }
                  } else if (is.null(recstr.t)) {
                    pred.measures['var.res.steps',h ,k ] <- Inf
                  }
                }
                step.chosen[, h] <- which(array(pred.measures['var.res.steps', , ],
                                                dim = c(n.steps, n.dims.loop) ) ==
                                          min(pred.measures['var.res.steps', , ], na.rm = TRUE),
                                          arr.ind = TRUE)[2:1]
                gapfill.results.step <- get(paste('gapfill.results.dim', step.chosen['dim', h], 
                                                  sep = ''))
                if (ratio.test.t != 1) {
                  slices.process     <- gapfill.results.step$slices.process
                  slices.constant    <- gapfill.results.step$slices.constant
                  values.constant    <- gapfill.results.step$values.constant 
                  slices.excluded    <- gapfill.results.step$slices.excluded
                  slices.without.gaps<- gapfill.results.step$slices.without.gaps
                  recstr.test        <- gapfill.results.step$reconstruction
                }
                if (print.status)
                  cat(paste(Sys.time(), ' : Chose dimension(s) ',
                            paste(dimensions[[ind]][[step.chosen['dim', h]]], collapse = ' and '), 
                            ' as first guess for next step.\n', sep = ''))
              } else {
                gapfill.results.step   <- gapfill.results
                gapfill.results.step$reconstruction[!is.na(recstr.test)] <-
                  recstr.test[!is.na(recstr.test)] 
              }
            } else {
              if (process.type == 'variances') {
                gapfill.results.step <- get(paste('gapfill.results.dim', step.chosen['dim', h], sep = ''))
              } else {
                gapfill.results.step <- gapfill.results
              }
            }        
          }
          
          ##determine first guess for next step
          if (!is.null(gapfill.results.step$reconstruction)) {
            results.reconstruction <- gapfill.results.step$reconstruction
            
            ## use first guess from other dimensions in case of too gappy series
            if (force.all.dims) {
              dim.other            <- setdiff(1:n.dims.loop, step.chosen['dim', h])
              if (length(dim.other) != 0) {
                results.dim.other  <- get(paste('gapfill.results.dim', dim.other, sep = ''))           
                if (!is.null(results.dim.other$reconstruction)) {
                  data.fill.otherdim     <- results.dim.other$reconstruction
                  data.fill.otherdim[!is.na(results.reconstruction)] <- NA

                  ## exclude not to be filled slices (oceans etc)
                  slices.fill.other <- !(gapfill.results.step$slices.process |
                                         gapfill.results.step$slices.too.gappy)               
                  if (sum(slices.fill.other) > 0) {
                    dim(slices.fill.other) <- dim(datacube)[gapfill.results.step$dims.cycle.id + 1]
                    ind.fill.other <- indexDatacube(datacube = datacube, logical.ind = slices.fill.other,
                                                   dims = gapfill.results.step$dims.cycle.id + 1)
                    data.fill.otherdim[ind.fill.other] <- NA
                  }
                  
                  n.data.fill.otherdim   <-  sum(!is.na(data.fill.otherdim))
                  if (n.data.fill.otherdim > 0) {
                    ratio <- round(n.data.fill.otherdim / sum(!is.na(results.reconstruction)) * 100, 2)
                    printStatus(paste('Including ',ratio, ' % values from dropped dimension ',
                                        paste(dimensions[[ind]][[dim.other]], collapse=','),
                                        sep = ''))
                    included.otherdim[h] <- TRUE
                    iterpath[dim(iterpath)[1], 'otherdim'] <- TRUE
                    results.reconstruction[is.na(results.reconstruction)] <-
                      data.fill.otherdim[is.na(results.reconstruction)]
                  }
                }              
              }
            }
            
            ##save first guess data
            if (drop.dim)
              results.reconstruction <- drop(results.reconstruction)
            file.name.guess.next   <- paste(sub('.nc$', '', file.name),
                                            '_first_guess_step_', formatC(h + 1, 2, flag = '0'),
                                            '_', process, '.nc', sep = '')
            file.copy(from = file.name, to = file.name.guess.next, overwrite = TRUE)          
            file.con.guess.next    <- open.nc(file.name.guess.next, write = TRUE)
            var.put.nc(file.con.guess.next, var.name, aperm(results.reconstruction, perm = aperm.reverse))
            close.nc(file.con.guess.next)
            step.use.frst.guess  <- min(c(h, step.chosen['step', max(which(!is.na(step.chosen['step',])))]))
            
          }
          ##TODO: add break criterium to get out of h loop
          ##      check what happens if gapfillSSA stops further iterations due to limiting groups of eigentriples
          
                                        # get iteration chosen information
          if (sum(!is.na(gapfill.results.step$iters.chosen)) > 0) {
            iter.chosen[1:2, h, process] <- apply(gapfill.results.step$iters.chosen, 2, mean, na.rm = TRUE)
            iter.chosen[3, h, process]   <- sum(is.na(gapfill.results.step$iters.chosen))
          }

          ##save process convergence information
          if (sum(!is.na(gapfill.results.step$process.converged)) > 0) {
            process.converged[process, h] <- sum(gapfill.results.step$process.converged, na.rm = TRUE) /
              sum(!is.na(gapfill.results.step$process.converged))
            process.converged.detail = c(process.converged.detail,
              list(array(gapfill.results.step$process.converged, dim = gapfill.results.step$dims.cycle.length)))
            name <- paste(process, '/', h, sep = '')
            names(process.converged.detail)[[length(process.converged.detail)]] <- name
            if (process.converged[process, h] < tresh.converged) {
              stop('More cells have not converged than allowed by tresh.converged!')
            }
          } 
        }
      }
      
      ##save results 
      .gapfillNcdfSaveResults(aperm.reverse = aperm.reverse, args.call.global = args.call.global,
                             datacube = datacube,
                             dims.cycle.id = gapfill.results.step$dims.cycle.id,
                             drop.dim = drop.dim,
                             file.name = file.name,                             
                             n.steps = n.steps,
                             ocean.mask = ocean.mask,
                             print.status = print.status, 
                             process.cells = process.cells,
                             reconstruction = results.reconstruction,
                             slices.without.gaps = gapfill.results.step$slices.without.gaps,                           
                             var.name = var.name,
                             var.names = var.names)
    }
    finished <- TRUE
    ## delete first guess files
    if (!keep.steps) {
      pattern.steps      <- paste(sub('[.]nc$', '_first_guess_step_', file.name))
      files.steps        <- list.files()[grepl(pattern.steps, list.files())] 
      file.remove(files.steps)
    }
    if (print.status)
      cat(paste(Sys.time(), ' : Gapfilling successfully finished. \n', sep = ''))
    close.nc(file.con.orig)
    if (process.type == 'variances') {
      out  <-list(pred.measures = pred.measures, step.chosen = step.chosen,
                  finished = finished, iterpath = iterpath, included.otherdim = included.otherdim,
                  SSAcallargs = args2SSA, iter.chosen = iter.chosen,
                  process.converged = process.converged, process.converged.detail =
                  process.converged.detail, seed = seed)
    } else {
      out  <- list(finished = finished, args2SSA = args2SSA, iterpath = iterpath,
                   iter.chosen = iter.chosen,
                   process.converged = process.converged, process.converged.detail =
                   process.converged.detail, seed = seed)
    }
    return(out)
  }, ex = function(){
    ## prerequisites: go to dir with ncdf file and specify file.name
    file.name        = 'scen_3_0.5_small.nc'
    max.cores        = 8
    calc.parallel    = TRUE
    
    ##example settings for traditional one dimensional "onlytime" setting and
    ## one step
    amnt.artgaps     <- list(list(c(0.05, 0.05)));
    amnt.iters       <- list(list(c(3, 10)));
    dimensions       <- list(list("time")); 
    M                <- list(list(12)); 
    n.comp           <- list(list(6)); 
    size.biggap      <- list(list(5)); 
    var.name         <- "auto"
    process.type     <- "stepwise"
#    .gapfillNcdf(file.name = file.name, dimensions = dimensions, amnt.iters = amnt.iters, 
#                amnt.iters.start = amnt.iters.start, amnt.artgaps = amnt.artgaps, 
#                size.biggap = size.biggap, n.comp = n.comp, tresh.fill = tresh.fill,
#                M = M, process.type = process.type)
    
    
    
    ##example settings for 3 steps, stepwise and mono dimensional
    dimensions       = list(list('time'), list('time'), list('time'))
    amnt.iters       = list(list(c(1,5)), list(c(2,5)), list(c(3,5)))
    amnt.iters.start = list(list(c(1,1)), list(c(2,1)), list(c(3,1)))
    amnt.artgaps     = list(list(c(0,0)), list(c(0,0)), list(c(0,0)))
    size.biggap      = list(list(0),      list(0),      list(0))
    n.comp           = list(list(6),      list(6),      list(6))
    M                = list(list(12),     list(12),     list(12))
    process.type     = 'stepwise'
#    gapfillNcdf(file.name = file.name, dimensions = dimensions, amnt.iters = amnt.iters, 
#                amnt.iters.start = amnt.iters.start, amnt.artgaps = amnt.artgaps, 
#                size.biggap = size.biggap, n.comp = n.comp, tresh.fill = tresh.fill,
#                M = M, process.type = process.type)
    
    ##example settings for 4 steps, stepwise and alternating between temporal and spatial
    dimensions       = list(list('time'), list(c('longitude','latitude')),
      list('time'), list(c('longitude','latitude')))
    amnt.iters       = list(list(c(1,5)), list(c(1,5)), list(c(2,5)), list(c(2,5)))
    amnt.iters.start = list(list(c(1,1)), list(c(1,1)), list(c(2,1)), list(c(2,1)))
    amnt.artgaps     = list(list(c(0,0)), list(c(0,0)), list(c(0,0)), list(c(0,0)))
    size.biggap      = list(list(0),      list(0), list(0),      list(0))
    n.comp           = list(list(15),     list(15), list(15),     list(15))
    M                = list(list(23),     list(c(20,20)), list(23), list(c(20,20)))
    process.type     = 'stepwise'
#    gapfillNcdf(file.name = file.name, dimensions = dimensions, 
#                amnt.iters = amnt.iters, amnt.iters.start = amnt.iters.start, 
#                amnt.artgaps = amnt.artgaps, size.biggap = size.biggap, n.comp = n.comp, 
#                tresh.fill = tresh.fill, M = M, process.type = process.type, max.cores = max.cores)
    
    ##example setting for process with alternating dimensions but variance criterium
    dimensions       = list(list('time', c('longitude','latitude')))
    n.comp           = list(list(5,      5))
    M                = list(list(10,     c(10, 10)))
    amnt.artgaps     = list(list(c(0,0), c(0,0)))
    size.biggap      = list(list(0,      0))
    process.type     = 'variances'
    max.steps        = 2
#    gapfillNcdf(file.name = file.name, dimensions = dimensions, n.comp = n.comp, 
#                tresh.fill = tresh.fill, max.steps = max.steps, M = M, 
#                process.type = process.type, amnt.artgaps = amnt.artgaps, 
#                size.biggap = size.biggap, max.cores = max.cores)
  })
globalVariables('i')


##################################### save results #############################
.gapfillNcdfSaveResults<- function(aperm.reverse, args.call.global, datacube, dims.cycle.id,
                                  drop.dim, file.name, n.steps,
                                  ocean.mask, print.status, process.cells,
                                  reconstruction, slices.without.gaps, var.name,
                                  var.names)
##title<< helper function for .gapfillNcdf
##details<< helper function for gapfillNcdf that saves the results ncdf file. 
##seealso<<
##\code{\link{gapfillNcdf}}  
{
  ##Prepare Results Ncdf File
  file.name.copy     <- paste(sub('[.]nc$','', file.name), '_gapfill.nc', sep = '')
  if (!file.exists(file.name.copy)) {
    copied             <- file.copy(from = file.name, to = file.name.copy, overwrite = TRUE)
    Sys.chmod(file.name.copy, mode = "0777")
    if (!copied)
      stop('Creating file for results failed!')
    file.name.copy     <- paste(sub('[.]nc$','', file.name), '_gapfill.nc', sep = '') 
    file.con.copy      <- open.nc(con = file.name.copy, write = TRUE)
    for (var.create in var.names) {
      var.def.nc(file.con.copy, paste(var.create, '_flag_orig', sep =''), 'NC_BYTE', 
                 var.inq.nc(file.con.copy, var.create)$dimids)
      att.put.nc(file.con.copy,  paste(var.create, '_flag_orig', sep =''),'long_name','NC_CHAR',
                 paste('flag indicating original values (1) and filled values (0) in ',
                       var.create, sep = ''))
      datacubeT                   <- transNcdfRotate(data.object = file.con.copy, var.name = var.name)$data
      data.flag                   <- array(NA, dim = dim(datacubeT))
      data.flag[is.na(datacubeT)]  <- 0
      data.flag[!is.na(datacubeT)] <- 1
      var.put.nc(file.con.copy, paste(var.name, '_flag_orig', sep =''),
                 aperm(data.flag,  perm = aperm.reverse))
      var.put.nc(file.con.copy, var.name, aperm(datacubeT,  perm = aperm.reverse))
    }  
  } else {
     file.con.copy      <- open.nc(con = file.name.copy, write = TRUE)     
  }
  
  ##prepare results
  dims.cycle.length               <- dim(datacube)[dims.cycle.id + 1]
  data.results.final              <- datacube
  data.results.final[is.na(data.results.final)] <- 
      reconstruction[is.na(data.results.final)]
  if (sum(slices.without.gaps) > 0 & process.cells == 'gappy') {
    dim(slices.without.gaps)      <- dims.cycle.length
    ind.array                     <- indexDatacube(datacube = datacube, 
                                                   logical.ind = slices.without.gaps, 
                                                   dims = dims.cycle.id + 1)
    data.results.final[ind.array] <- datacube[ind.array]
  }
  if (length(ocean.mask) > 0) {
    ind.array                     <- indexDatacube(datacube = datacube, 
                                                  logical.ind = ocean.mask, 
                                                  dims = c(1,2))
    data.results.final[ind.array] <- NA
  }

  #save results
  if (print.status)
    cat(paste(Sys.time(), ' : Writing results to file. \n', sep = ''))
  if (drop.dim)
    data.results.final <- drop(data.results.final)
  var.put.nc(file.con.copy, var.name, aperm(data.results.final,  perm = aperm.reverse))
  sync.nc(file.con.copy)

  #add attributes with process information to ncdf files
  if (var.name == var.names[length(var.names)]) {
    string.args  <- args.call.global
    att.put.nc(file.con.copy, 'NC_GLOBAL', 'Filled_by', 'NC_CHAR', getSysinfo())
    att.put.nc(file.con.copy, 'NC_GLOBAL', 'Filled_on', 'NC_CHAR', as.character(Sys.time()))
    att.put.nc(file.con.copy, 'NC_GLOBAL', 'Filling_settings', 'NC_CHAR', string.args)
    hist.string.append <- paste('Gaps filled on ', as.character(Sys.time()), ' by ',
                                Sys.info()['user'], sep = '')
    if (is.element('history', infoNcdfAtts(file.con.copy, 'NC_GLOBAL')[, 'name'])) {
      hist.string    <- paste(att.get.nc(file.con.copy, 'NC_GLOBAL', 'history'), 
                              '; ', hist.string.append)
      att.put.nc(file.con.copy, 'NC_GLOBAL', 'history', 'NC_CHAR', hist.string)
    } else {
      att.put.nc(file.con.copy, 'NC_GLOBAL', 'history', 'NC_CHAR', hist.string.append)
    }
  }
  close.nc(file.con.copy)
}


##################################   gapfill data cube #########################
.gapfillNcdfDatacube <- function(args.call.SSA = list(), calc.parallel = TRUE,
                                datacube, datapts.n,  debugging, dims.info, dims.cycle,
                                dims.cycle.id,dims.process,  dims.process.id,
                                file.name, first.guess = 'mean', g, h, iters.n,
                                l, max.cores = 8,  MSSA, MSSA.blocksize,
                                MSSA.blck.trsh = MSSA.blck.trsh, ocean.mask = c(),
                                print.status, process.cells = c('gappy','all')[1], 
                                ratio.const = ratio.const, ratio.test.t, reproducible,
                                slices.constant = c(), slices.excluded = c(),
                                slices.process = c(), slices.without.gaps= c(),
                                tresh.const = tresh.const, tresh.fill.dc =  .1,
                                values.constant = c())
##title<< helper function for gapfillNcdf
##details<< helper function for gapfillNcdf that handles the main datacube transformations. 
##seealso<<
##\code{\link{gapfillNcdf}}    
{
  slices.n            <- prod(dim(datacube)[dims.cycle.id + 1])
  dims.process.length <- dim(datacube)[dims.process.id + 1]
  dims.cycle.length   <- dim(datacube)[dims.cycle.id + 1]
  
  #identify valid grids
  args.identify <- list(dims.cycle = dims.cycle, dims.cycle.id = dims.cycle.id,
                        dims.process.id = dims.process.id, datacube = datacube,
                        MSSA = MSSA, dims.process =  dims.process, 
                        process.cells = c('gappy','all')[1], first.guess = first.guess,
                        ocean.mask = ocean.mask, print.status = print.status, 
		                    slices.n = slices.n, dims.process.length = dims.process.length,
                        tresh.fill.dc = tresh.fill.dc, ratio.test.t = ratio.test.t,
                        g = g, ratio.const = ratio.const, algorithm = 'Gapfill', 
                        tresh.const = tresh.const, args.call.SSA = args.call.SSA)
  if (g == 2)
    args.identify <- c(args.identify, list(slices.excluded = slices.excluded,
                                           values.constant = values.constant,
                                           slices.constant = slices.constant,
                                           slices.process = slices.process))
  results.identify      <- do.call('.identifyValidCellsSSA', args.identify)
  iters.n = results.identify$iters.n
  slices.process = results.identify$slices.process
  values.constant = results.identify$values.constant
  slices.constant = results.identify$slices.constant 
  slices.without.gaps = results.identify$slices.without.gaps
  slices.excluded = results.identify$slices.excluded
  slices.too.gappy = results.identify$slices.too.gappy
  args.call.SSA$tresh.min.length =  results.identify$tresh.min.length
  
  #create iterator
  if (sum(slices.process) == 0) {
    data.results.finished <- array(NA, dim(datacube))
    data.variances        <- NULL
    iters.chosen          <- c(NA, NA)
    process.converged     <- rep(NA, prod(dims.cycle.length))
  } else {
    results.crtitercube  <- do.call(.gapfillNcdfCreateItercube, 
                                    list(datacube = datacube, 
                                         slices.process = slices.process,
                                         iters.n = iters.n, slices.n = slices.n, 
                                         dims.cycle.length = dims.cycle.length, 
                                         dims.cycle.id = dims.cycle.id, 
                                         max.cores = max.cores, MSSA = MSSA, 
                                         MSSA.blocksize = MSSA.blocksize,
                                         MSSA.blck.trsh = MSSA.blck.trsh))
    iter.gridind = results.crtitercube$iter.gridind
    ind.process.cube = results.crtitercube$ind.process.cube
    max.cores = results.crtitercube$max.cores
    index.MSSAseries = results.crtitercube$index.MSSAserie

    
    #perform (parallelized) calculation
    if (print.status)
      cat(paste(Sys.time(), ' : Starting calculation: Filling ', sum(slices.process),
              ' time series/grids of size ', datapts.n, '. \n', sep = ''))

    if (calc.parallel) {
      results.parallel = foreach(i = 1:max.cores                       
        , .combine = rbindMod
        , .multicombine = TRUE, .packages = 'spectral.methods') %dopar% gapfillNcdfCoreprocess(
                                  iter.nr = i, datacube = datacube,
                                  dims.process.id = dims.process.id, datapts.n = datapts.n, args.call.SSA = args.call.SSA,
                                  iter.gridind = iter.gridind, ind.process.cube = ind.process.cube, first.guess = first.guess,
                                  print.status = print.status, iters.n = iters.n, dims.cycle.length = dims.cycle.length, 
                                  dims.cycle.id = dims.cycle.id,  dims.process.length =  dims.process.length, MSSA = MSSA, 
                                  file.name = file.name, reproducible = reproducible)
    } else {
       results.parallel = foreach(i = 1:1
         , .combine =  rbindMod
         , .multicombine = TRUE, packages = 'spectral.methods') %do% gapfillNcdfCoreprocess(iter.nr = i, datacube = datacube,
                                   dims.process.id = dims.process.id, datapts.n = datapts.n, args.call.SSA = args.call.SSA,
                                   iter.gridind = iter.gridind, ind.process.cube = ind.process.cube, first.guess = first.guess,
                                   print.status = print.status, iters.n = iters.n, dims.cycle.length = dims.cycle.length, 
                                   dims.cycle.id = dims.cycle.id,  dims.process.length =  dims.process.length, MSSA = MSSA, 
                                   file.name = file.name, reproducible = reproducible)            
     }
    
    data.results.valid.cells <- results.parallel$reconstruction
    data.variances           <- results.parallel$variances
    iters.chosen             <- results.parallel$iters.chosen
    process.converged        <- array(NA, dim = dims.cycle.length)
    ## if (interactive()) {
    ##   #browser()
    ## } else {
    ##   path.save = '/Net/Groups/BGI/people/jbuttlar/Scratch'
    ##   file.save    <- paste('debug_', file.name, '_step_1.RData',  sep = '')
    ##   printStatus(paste('Saving debugging file to ', file.path(path.save, file.save)))
    ##   save(list = ls(),  file = file.path(path.save, file.save))
    ## }
    process.converged[slices.process]  <- results.parallel$process.converged

    #fill all values to results array
    if (print.status)
      cat(paste(Sys.time(), ' : Transposing results. \n', sep = ''))
    data.results.all.cells                    <- array(NA, dim = c(slices.n, datapts.n))
    if (sum(slices.process) > 0)
      data.results.all.cells[index.MSSAseries, ]<- data.results.valid.cells
    data.results.all.cells[slices.constant, ] <- rep(values.constant[slices.constant], times = datapts.n)
    
    #reshape results array to match original data cube
    data.results.rshp          <- array(data.results.all.cells, dim = c(dims.cycle.length, dims.process.length))
    dims.order.results         <- c(dims.cycle, dims.process)
    if (dim(datacube)[1] == 1 && length(dim(datacube)) == 2) {                      ## if only one series is existent
      perm.array = c(1,2)
    } else {
      perm.array                 <- match(dims.info[, 'name'], dims.order.results)
    }
    data.results.finished      <- aperm(data.results.rshp, perm.array)
  }

  return(list(reconstruction = data.results.finished, data.variances = data.variances,
              slices.without.gaps = slices.without.gaps, slices.process = slices.process, 
              max.cores = max.cores, slices.process = slices.process,
              slices.too.gappy = slices.too.gappy, sices.constant = slices.constant,
              values.constant = values.constant, slices.excluded = slices.excluded,
              iters.chosen = iters.chosen, process.converged = process.converged))
}



###############################  create iteration datacube   ###################


.gapfillNcdfCreateItercube  <- function(datacube, iters.n, dims.cycle.length, 
    dims.cycle.id, slices.process, max.cores, slices.n, MSSA, MSSA.blocksize,
    MSSA.blck.trsh)
##title<< helper function for gapfillNcdf
##details<< helper function for gapfillNcdf that creates the index array used
##          in the foreach iterations to extract data.     
##seealso<<
##\code{\link{gapfillNcdf}}      
{
  ##TODO
  #make indices independent from dimension order
  index.MSSAseries   <- integer()
  index.MSSAnr       <- array(1:prod(dims.cycle.length), dim = dims.cycle.length)
  if (MSSA) {
    slices.remove    <- array(!slices.process, dim = dims.cycle.length)
    iters.n          <- prod(ceiling(dim(datacube)[dims.cycle.id + 1] / MSSA.blocksize))
    ind.process.cube <- array(FALSE, dim = c(iters.n, dims.cycle.length))
    row.block        <- 1
    col.block        <- 1
    for (p in 1:iters.n) {      
      ind.row        <- row.block:(min(c(dim(datacube)[1], (row.block + MSSA.blocksize - 1))))
      ind.col        <- col.block:(min(c(dim(datacube)[2], (col.block + MSSA.blocksize - 1))))
      mn.blck.fll    <- MSSA.blck.trsh * prod(length(ind.row), length(ind.col),
                                              dim(datacube)[-(dims.cycle.id + 1)])
      if (sum(!is.na(datacube[ind.row,ind.col,])) > mn.blck.fll) 
        ind.process.cube[p, ind.row, ind.col] <- TRUE
      #create vector with indices for the time series computed
      ind.add          <- index.MSSAnr[ind.process.cube[p,,]][!slices.remove[ind.process.cube[p,,]]]
      index.MSSAseries <- c(index.MSSAseries, ind.add)
      #determine indices for next iter
      row.block        <- row.block + MSSA.blocksize      
      if (row.block > dim(datacube)[1]) {
        row.block      <- 1
        col.block      <- col.block + MSSA.blocksize
      }         
    }
    
    ##remove ocean cells
    if (sum(!slices.process) > 0) {
      ind.process.cube[indexDatacube(datacube = ind.process.cube, 
                                    logical.ind = slices.remove, dims = c(2, 3))] <- FALSE
    }
    rows.remove <- apply(ind.process.cube, 1, function(x){sum(x) == 0})
    if (sum(rows.remove) > 0) {
      ind.process.cube <- ind.process.cube[!rows.remove, , ]
      iters.n <- iters.n - sum(rows.remove)
    } 
            
  } else if (!MSSA) {
    args.expand.grid     <- alist()
    for (d in 1:length(dims.cycle.id))
      args.expand.grid[[d]] <- 1:dim(datacube)[dims.cycle.id[d] + 1]    
    iter.grid.all <- cbind(1:slices.n, as.matrix(do.call("expand.grid",
                                                         args.expand.grid)))
    ind.process.cube     <- array(iter.grid.all[slices.process, ],
                                  dim = c(sum(slices.process), length(dims.cycle.id) + 1))
    index.MSSAseries     <- (1:slices.n)[slices.process]
  }
  max.cores              <- min(c(iters.n, max.cores))
  iters.per.cyc          <- rep(floor(iters.n / max.cores), times = max.cores)
  if (!(iters.n %% max.cores) == 0)
    iters.per.cyc[1:((iters.n %% max.cores))] <- floor(iters.n / max.cores) + 1
  iter.gridind           <- matrix(NA, ncol = 2, nrow = max.cores)
  colnames(iter.gridind) <- c('first.iter', 'last.iter')
  if (max.cores > 1)  {
    iter.gridind[, 1]    <- (cumsum(iters.per.cyc) + 1) - iters.per.cyc
    iter.gridind[, 2]    <- cumsum(iters.per.cyc)
  } else {
    iter.gridind[]       <- c(1, iters.n)
  }
  return(list(iter.gridind = iter.gridind , ind.process.cube = ind.process.cube,
              max.cores = max.cores, index.MSSAseries = index.MSSAseries))
}


##################  combine data from foreach iteration ########################
rbindMod <- function(...) 
##title<< helper function for gapfillNcdf
##details<< helper function for gapfillNcdf that combines the foreach output
##          in a convenient way.
{
  
  cat(paste(Sys.time(), ' : Assembling data from parallelized computations.\n', 
            sep=''))
  assign('dummy', list(...))


  path.save = '/Net/Groups/BGI/people/jbuttlar/Scratch'
  file.save    <- paste('debug_combine_', as.numeric(Sys.time()), '_', sample(1:100, 1), '.RData',  sep = '')
  printStatus(paste('Saving debugging file to ', file.path(path.save, file.save)))
  save(dummy,  file = file.path(path.save, file.save))
  


  

  vars.amnt <- dim(dummy[[1]][['variances']])[2]
  cube.cols <- sum(sapply(dummy,function(x)dim(x[['reconstruction']])[1]))
  reconstruction <- matrix(unlist(lapply(dummy, function(x)as.vector(t(x[[1]])))),
                           nrow = cube.cols, byrow = TRUE)
  variances      <- matrix(unlist(lapply(dummy, function(x)as.vector(t(x[[2]])))),
                           ncol = vars.amnt, byrow = TRUE)
  process.converged <- unlist(lapply(dummy, function(x)as.vector(t(x[['process.converged']]))))
  iters.chosen   <- matrix(unlist(lapply(dummy, function(x)as.vector(t(x[['iters.chosen']])))),
                           ncol = 2, byrow = TRUE)
  
  return(list(reconstruction = reconstruction, variances = variances, 
              process.converged = process.converged, iters.chosen = iters.chosen))    
}


########################## gapfill function for single core ####################
gapfillNcdfCoreprocess <- function(args.call.SSA, datacube, datapts.n, dims.cycle.id,
                                   dims.cycle.length, dims.process.id,
                                   dims.process.length, file.name, first.guess,
                                   ind.process.cube, iter.gridind, iter.nr, iters.n,
                                   MSSA, print.status, reproducible)
##title<< helper function for gapfillNcdf
##details<< helper function for gapfillNcdf performs each individual series/grid 
##          extracion and handing it over to gapfillSSA.
{
   ## path.save = '/Net/Groups/BGI/people/jbuttlar/Scratch'
   ## file.save <- file.path(path.save, paste('debug_', file.name, '_SSA_internal_core_single_',
   ##                                         args.call.SSA$iterinf$process, '_',
   ##                                         args.call.SSA$iterinf$h,'_', args.call.SSA$iterinf$l, '_', 
   ##                                          iter.nr, '.RData',  sep = ''))
   ## printStatus(paste('Saving debugging file to ', file.save))
   ## save(list = ls(environment()), file = file.save)

  
  #determine iteration parameters
  iter.ind          <- iter.gridind[iter.nr, ]
  datapts.n         <- prod(dim(datacube)[dims.process.id + 1])
  n.series.steps    <- numeric()
  if (MSSA){
    ind.t           <- rep(FALSE, times = dim(ind.process.cube)[1])
    ind.t[iter.ind[1]:iter.ind[2]]       <- TRUE
    n.series.core   <- sum(ind.process.cube[indexDatacube(datacube = ind.process.cube, 
                                                         logical.ind = ind.t, dims = 1)])
  } else {
    n.series.core   <- (diff(iter.ind) + 1)
  }

  # define results and diagnostics arrays
  data.results.iter <- array(NA, dim = c(n.series.core , datapts.n))
  variances         <- array(NA, dim = c(diff(iter.ind) + 1, args.call.SSA[['n.comp']]))
  process.converged <- array(NA, dim = c(diff(iter.ind) + 1))
  iters.chosen      <- array(NA, dim = c(diff(iter.ind) + 1, 2))
  
  for (n in 1:(diff(iter.ind) + 1)) {
    ind.total       <- rep(FALSE, iters.n)
    ind.total[(iter.ind[1] : iter.ind[2])[n]] <- TRUE 
    data.results.iter.t = try({

      ## determine arguments transferred to SSA process
      args.call.t             <- args.call.SSA
      args.call.t$iterinf   =  NULL
      dims.extr.data          <- dims.process.length
      aperm.extr.data         <- 1:(length(dims.process.id) + 1) 
      if (MSSA) {
        ind.act.cube          <- array(ind.process.cube[indexDatacube(datacube = ind.process.cube, 
                                                                     logical.ind = ind.total, 
                                                                     dims = 1)],
                                       dim = dims.cycle.length)
        n.series.steps[n]     <- sum(ind.act.cube)
        dims.extr.data        <- c(dims.extr.data, n.series.steps[n])
        aperm.extr.data       <- c(2,1)
        perm.before           <- 1:2
        ind.extr              <- indexDatacube(datacube = datacube, 
                                              logical.ind = ind.act.cube, 
                                              dims = dims.cycle.id + 1)
      } else {
        perm.before           <- order(dims.process.id)
        ind.matrix            <- array(FALSE, dims.cycle.length)
        ind.matrix.list       <- matrix(ind.process.cube[(iter.ind[1] : iter.ind[2])[n], -1], byrow=TRUE,
                                        ncol = length(dims.cycle.length))
        ind.matrix[ind.matrix.list] <- TRUE
        ind.extr              <- indexDatacube(datacube = datacube, 
                                              logical.ind = ind.matrix, 
                                              dims = dims.cycle.id + 1)
        n.series.steps[n]     <- 1
      }     
      series.noperm                  <- array(datacube[ind.extr], dim =  dims.extr.data)
      args.call.t[['series']]        <- aperm(series.noperm, perm = perm.before)
      if (!(class(first.guess) == 'character' && first.guess == 'mean')) {
        fg.noperm                    <- array(first.guess[ind.extr], dim =  dims.extr.data)
        args.call.t[['first.guess']] <- aperm(fg.noperm, perm = perm.before)         
      }

      ## run SSA
      series.filled       <- do.call(gapfillSSA, args.call.t)
      
      ## transpose and extract SSA results
      rcstr.local         <- aperm(array(series.filled$reconstr,
                                         dim = c(dims.process.length, n.series.steps[n])),
                                   aperm.extr.data)
      ind.results <- (1 : n.series.steps[n]) + (((n>1) * 
                                                 (sum( n.series.steps[1 : max(c(n - 1, 1))]))))  
      data.results.iter[ind.results, ]  <- array(rcstr.local, dim = c(n.series.steps[n], datapts.n))
      variances[n, 1:length(as.vector(series.filled$variances))] <- as.vector(series.filled$variances)
      process.converged[n]              <- series.filled$process.converged
      iters.chosen[n, ]                 <- series.filled$iter.chosen
      'completed'      
    })
   # browser()
    
    
    ## save workspace in case of error for debugging
    if (class(data.results.iter.t) == 'try-error') {
      print(paste('Error occoured at iteration ', iter.nr, ' and loop ', n, '!', sep = ''))
      error.from.calc                 <- data.results.iter.t
      data.results.iter.t             <- matrix(Inf, ncol = datapts.n, nrow = 1)
      system.info                     <- sessionInfo()    
      path.file                       <- file.path('/', 'Net', 'Groups', 'BGI', 'tmp', 
                                                   'jbuttlar', 'Cluster_jobs_debugging', getwd()) 
      file.name.t                     <- file.path(path.file, paste('workspace_error_', file.name, '_',
                                                                    iter.nr, '_', n, sep = ''))
      dump.frames(to.file = TRUE, dumpto = file.name.t)
      print(paste('Saved workspace to file ', file.name.t, '.rda', sep = ''))
      stop()
    }

    ## print loop progress information
    if (iter.nr == 1 &&( diff(iter.ind) < 20  || 
          (n %% (ceiling((diff(iter.ind)) / 20)) == 0)))
      if (print.status)
        cat(paste(Sys.time(), ' : Finished ~',
                  round(n / (diff(iter.ind) + 1) * 100), '%. \n', sep=''))
  }
  


  



  
  ## return results
  return(list(reconstruction = data.results.iter, variances = variances, 
              process.converged = process.converged, iters.chosen = iters.chosen))
}

#####################       # check input to functions      ####################

.checkInputNcdfSSA <- function(
  SSAprocess ##<< character string: name of the function for which to check the
             ##   arguments. One of 'Gapfill' or 'Decompose'.
  , ...)
  ##title<< check input parameters for the SSA  routines.
  ##details<< helper function for Ncdf SSA (Decompose, Gapfill) routines that
  ##          checks the consistency of the input parameters.
  ##seealso<<
  ##\code{\link{gapfillNcdf}}, \code{\link{decomposeNcdf}}
{    
  ##TODO
  # include test for MSSA and windowlength=1
  
  ##TODO
  #check input file
  ##TODO add checks
  ## - single dimension variances and thresh.fill.first, tresh.fill
  ## - facilitate old school SSA via gaps.cv, max.steps = 1, amnt.artgaps !=0
  ##   and length(dimensions) == 1
  ##TODO check intercorelation between ratio.test and gaps.cv
  args <- list(...)
      
  if (is.null(args$file.name) )
    stop('file.name needs to be supplied!')
  
  if (!file.exists(args$file.name))
    stop('Input ncdf file not found. Check name!')
  if (SSAprocess == 'Decompose') {
     dims.check  <- 'time'
  } else {
     dims.check <- unique(unlist(args$dimensions))
  }  
  check.passed <- checkNcdfFile(file.name = args$file.name, 
                                 dims = dims.check, 
                                 type = 'relaxed')
  if (!check.passed)
    stop('NCDF file not consistent with CF ncdf conventions!')
  file.con.orig <- open.nc(args$file.name)
  if (args$var.names[1] != 'auto')
    for (var.name in args$var.names)
      if (!is.element(var.name, infoNcdfVars(file.con.orig)$name))
        stop(paste('Variable name ', var.name, 'does not exist in ncdf file!', sep = ''))
  
  if(sum(is.na(match(unique(unlist(args$dimensions)), c('longitude', 'latitude', 'time')))) > 0)
    stop('Every dimensions value has to be one of time, longitude, latitude!')
  dims.not.exist <- is.na(match(unlist(args$dimensions), infoNcdfDims(file.con.orig)[, 'name']))
  if(sum(dims.not.exist) > 0)
    stop(paste('Dimension(s) ', paste(unlist(args$dimensions)[dims.not.exist], collapse = ', '),
            'not existent in ncdf file!', sep = ''))

  if (SSAprocess == 'Decompose' ) {
    n.bands                       <- length(unlist(args$borders.wl))-length(args$borders.wl)
    if(file.info(args$file.name)$size / 1024^3 * n.bands >  2)
      stop('Target file size may exceed 2GB. Reduce n.bands or compress or split input ncdf file!')
    if (is.null(args$borders.wl))
      stop('Argument borders.wl needs to be supplied!')
    if (!class(args$borders.wl) == 'list')
      stop('Wrong class for borders.wl! Needs to be a list!')
    args.return    <- list(file.con.orig = file.con.orig)    
  } else if (SSAprocess == 'Gapfill') {
    ## check first guess file
    name.first.guess.std <- paste(sub('.nc', '', args$file.name),'_first_guess_step_1.nc', sep = '')
    if (!(args$first.guess == 'mean')) {
      if (!(file.exists(args$first.guess))) {
        stop('Specified first.guess file not existent!')
      } else if (!(args$first.guess == name.first.guess.std)) {
        stop('Name for supplied first guess file does not fit the standardized scheme!')
      }
      check.passed = checkNcdfFile(file.name = args$first.guess, 
        dims = infoNcdfDims(file.con.orig)[, 'name'], 
        type = "relaxed")
      if (!check.passed)
        stop('First guess NCDF file not consistent with CF ncdf conventions!')
    }
    
    ##transfer and check ocean mask
    lengths.dim.nontime <- infoNcdfDims(file.con.orig)['length'][!(infoNcdfDims(file.con.orig)['name'] == 'time')]
    
    if (!is.null(args$ocean.mask) && class(args$ocean.mask) == 'character') {
      if (!(file.exists(args$ocean.mask)))
        stop('File for ocean mask not existent!')
      check.passed <- checkNcdfFile(file.name = args$ocean.mask, 
                                      dims=c('longitude', 'latitude'), type = 'relaxed')
      if (!check.passed)
        stop('ocean mask NCDF file not consistent with CF ncdf conventions!')
      con.ocean   <- open.nc(args$ocean.mask)
      var.names.t <- infoNcdfVars(con.ocean)[, 'name']
      var.name.om <- var.names.t[is.na(match(var.names.t, c('time', 'longitude', 'latitude')))]
      if (length(var.name.om) > 1)
        stop('More then one variable existent in ocean mask!')
      for (par.check in c('longitude', 'latitude')) {
        if(!identical(var.get.nc(file.con.orig, par.check), 
                      var.get.nc(con.ocean, par.check)))
          stop(paste(par.check,' values need to be identical in ocean.mask and',
                     'input ncdf file!'))
      }
      ocean.cells <- var.get.nc(con.ocean, var.name.om)
      args$ocean.mask  <- array(FALSE, dim = dim(ocean.cells))
      args$ocean.mask[ocean.cells == 1 ] <- TRUE
      data.orig      <- var.get.nc(file.con.orig, readNcdfVarName(file.con.orig))
      oceancells.data<- sum(!is.na(data.orig[indexDatacube(data.orig, args$ocean.mask, c(1,2))]))
      if (oceancells.data > 0)
        stop('Some ocean cells seem to contain data. Is the ocean.mask file correctly set up?')
      close.nc(con.ocean)
    }
    if (length(args$ocean.mask) > 0 && !(dim(args$ocean.mask) == lengths.dim.nontime))
    stop(paste('The ocean mask has to have identical dimensions as the spatial',
               ' dimensions in the ncdf file!', sep = ''))
    
    ## check consitency of input   
    n.steps           <- length(args$amnt.artgaps)
    if (!is.element(args$process.cells, c('gappy', 'all')))
      stop('process.cells has to be one of \'all\' or \'gappy\'!')
    if (!is.element(args$process.type, c('stepwise', 'variances')))
      stop('process.type has to be one of \'stepwise\' or \'variances\'!')
    if (args$max.cores > 1 & !args$calc.parallel)
      stop('More than one core can only be used if calc.parallel==TRUE!')
    if (!inherits(args$tresh.fill, 'numeric') || length(args$tresh.fill) != 1)
      stop('Supply argument tresh.fill as object of class numeric and length = 1!')
    if (args$tresh.fill < 0 | args$tresh.fill > 1)
      stop('All values in tresh.fill need to be between 0 and 1!')

    dims.order <- c('latitude', 'longitude', 'time')
    for (i in 1:length(args$dimensions))
      for (j in 1:length(args$dimensions[[i]])) {
        dims.check <-  args$dimensions[[i]][[j]]
        if (length(dims.check) > 1)
          if (sum(diff(match(args$dimensions[[i]][[j]],  dims.order)) < 0) > 0)
            stop(paste('Datacube is internally rotated/transposed to the dimension',
                       'oder: lat/lon/time. Please supply all dim pairs in dimensions',
                       'in this order (i.e. lat/lon,  lon/time or lat/time.).'))
      }
    if (args$process.type == 'variances') {
      args.check <- c('dimensions', 'n.comp', 'M', 'amnt.artgaps', 'size.biggap')
      n.dims     <- length(args$dimensions[[1]])
      for (m in 1:length(args.check)) {
        if (length(args[[args.check[m]]][[1]]) != n.dims)
          stop(paste(args.check[m],' needs to be of the same length as dimensions',
                     '[[1]] for process.type == \'variances\'!', sep = ''))
        if(!(class(args[[args.check[m]]])) == 'list')
          stop(paste('Argument ',args.check[m], ' is not a list!'))
        if (length(args[[args.check[m]]]) != 1)
          stop(paste(args.check[m],' needs to be of length 1 for process.type ==',
                   ' \'variances\' !', sep = ''))
      }
    }  else if (args$process.type == 'stepwise') {
      args.list = c('amnt.artgaps', 'size.biggap', 'n.comp', 'M', 'pad.series', 
         'amnt.iters', 'amnt.iters.start', 'MSSA')
      for (n in 1:length(args.list)) {
        if(class(args[[args.list[n]]]) != 'list')
          stop(paste('Argument ', args.list[n], ' is not a list!'))
        if (length(args[[args.list[n]]]) != length(args$dimensions))
          stop(paste('Argument ', args.list[n], ' has to be a list of the same ',
                     'length as dimensions for process.type == \'stepwise\'!'))
        ##      if (! all(sapply(get(args.list[n]), length) == 1 ))
        ##        stop(paste('Argument ', args.list[n], '[[1..n]] can only be a list of ',
        ##                   'length one for process.type == \'stepwise\'!', sep = ''))
      }
      for (i in 1:length(args$MSSA)) {
        for (j in 1:length(args$MSSA[i]))
          if (args$MSSA[[i]][[j]] & length(args$M[[i]][[j]]) != 2)
            stop('If MSSA ought to be computed, all corresponding Ms need to be of length 2.')
      } 
      for (o in 1:length(args$dimensions))
        if (sapply(args$dimensions[o] , function(x) length(x)) > 2)
          stop('Settings imply SSA with more than two dimensions which is not implemented.')
      step.wrong <- (unlist(args$tresh.fill) == 0 & (!(unlist(args$size.biggap) == 0) | 
                             !(all(unlist(args$amnt.artgaps) == 0))))
    }
    close.nc(file.con.orig)
    args.return <- list(ocean.mask = args$ocean.mask)
    
  }
  return(invisible(args.return))  
}

########################       identify valid cells         ####################

.identifyValidCellsSSA <- function(
##title<< helper function for Ncdf SSA algorithms 
    dims.cycle.id, dims.process.id, datacube, ratio.const, tresh.const , print.status,
    slices.n, algorithm,  g = c() ,process.cells = c('gappy','all')[1], dims.cycle = c(), 
    args.call.SSA = list(), tresh.fill.dc = 0, ratio.test.t =1, first.guess = 'mean', 
    ocean.mask = c(), dims.process = c(), dims.process.length = 0, 
    slices.without.gaps = rep(FALSE, slices.n), slices.too.gappy = rep(FALSE, slices.n),
    slices.constant = rep(FALSE, slices.n), slices.process = rep(TRUE, slices.n),
    slices.ocean = rep(FALSE, slices.n), values.constant = integer(length = slices.n),
    slices.excluded = rep(FALSE, slices.n), MSSA = FALSE)   
##description<<
##  helper function for gapfillNcdf and decomposeNcdf that identifies the grid cells to process. 
##seealso<<
##\code{\link{gapfillNcdf}}, \code{\link{gapfillNcdf}}    
{  
  ## TODO
  ## -possibility to identify gap less MSSA blocks
  ## -include possibility to infer slices.continuous max border

  if (print.status)
    cat(paste(Sys.time(), ' : Identifying valid cells ...\n', sep=''))
  treshhld.gappy <-  1

  ## get amount of missing values
  if (algorithm == 'Gapfill') {
    add.id = 1  
  } else {
    add.id = 0
  }  
  getMissingRatio = function(x) {
    sum(is.na(x)) / prod(dim(datacube)[dims.process.id + add.id])
  } 
  amnt.na      <- apply(datacube, MARGIN = dims.cycle.id + add.id , getMissingRatio)
  slices.empty <- as.vector(amnt.na == 1)

  ## slices checks specific to decomposition
  if (algorithm == 'Decompose') { 
    slices.process              <- as.vector(amnt.na == 0)
    slices.too.gappy            <- !slices.empty & !slices.process
    if (sum(slices.too.gappy) > 0) {
       slices.continuous         <- apply(datacube, MARGIN = dims.cycle.id + add.id , isSeriesContinuous)
       slices.continuous[(1-amnt.na[slices.continuous]) * dim(datacube)[dims.process.id] < 20] <- FALSE
       slices.too.gappy[slices.continuous] <- FALSE
    }
    if (sum(slices.too.gappy) > 0)
      cat(paste(Sys.time(), ' : ', sum(slices.too.gappy),' series with gaps were found. ',
                            'Spectral decomp. for these is not possible!\n',sep=''))
 
    ## slices checks specific to gap filling 
  } else if (algorithm == 'Gapfill') {
    if (!MSSA) {
      if (g == 1) {      
        ## modifications in case of given ocean mask
        if ((length(ocean.mask) > 0 ) & (sum(!is.na(match(c('longitude','latitude'), 
                          dims.process))) == 2)) {
          amnt.na <- 1- apply(datacube, MARGIN = dims.cycle.id + 1 ,
              function(x) sum(!is.na(x[as.vector(!ocean.mask)])) / sum(!ocean.mask)   )
        }
        if (length(ocean.mask) > 0 & sum(!is.na(match(c('longitude','latitude'), 
                    dims.cycle))) == 2) 
          slices.ocean           <- as.vector(ocean.mask)
        
        ## exclude non gappy slices
        if (process.cells == 'gappy') 
          slices.without.gaps       <- as.vector((amnt.na == 0))
        
        ## exclude too gappy slices
        gaps     <- args.call.SSA$amnt.artgaps
        size.bg  <- args.call.SSA$size.biggap^(length(dims.process.length))
        dtpts    <- prod(dims.process.length)
        if (sum(gaps != 0) > 0 & tresh.fill.dc > 0) {
          n.biggaps   <- max(c(floor(dtpts * gaps[1] /size.bg), 1))  
          if (gaps[1] > 0) {
            n.smallgaps <- n.biggaps * size.bg * gaps[1] / gaps[2]
          } else {
            n.smallgaps <- dtpts * gaps[2]
          }
          treshhld.gappy <- 1 - ((n.smallgaps + n.biggaps * size.bg) / dtpts ) - tresh.fill.dc
        } else {
          treshhld.gappy <- 1 - (tresh.fill.dc)
        }
        slices.too.gappy <- as.vector(amnt.na >= treshhld.gappy)   
        slices.too.gappy[slices.ocean] <- FALSE             
      } else if (g != 1) {
        slices.process <- slices.excluded
        iters.n <- sum(slices.process)
        return(list(iters.n = iters.n, slices.process = slices.process, 
                values.constant = values.constant, slices.constant = slices.constant, 
                slices.without.gaps = slices.without.gaps, 
                slices.excluded = slices.excluded,
                slices.too.gappy = slices.too.gappy)) 
      }
    } else if (MSSA) {
      if (length(ocean.mask) > 0) {
        slices.ocean                 <- as.vector(ocean.mask)
        slices.process[slices.ocean] <- FALSE
        iters.n <- sum(slices.process)
        return(list(iters.n = iters.n, slices.process = slices.process, 
                    values.constant = values.constant, slices.constant = slices.constant, 
                    slices.without.gaps = slices.without.gaps, 
                    slices.excluded = slices.excluded,
                    slices.too.gappy = slices.too.gappy)) 
      }
    }   
  }
  
  # identify constant slices
  slices.constant    <- as.vector(apply(X = datacube, MARGIN= dims.cycle.id + add.id,
                                        FUN = isSeriesConstant, ratio.const = ratio.const,
                                        tresh.const = tresh.const))
  slices.constant[slices.too.gappy | slices.empty | slices.ocean | slices.without.gaps] <- FALSE
  values.constant    <-  as.vector(apply(datacube, MARGIN = dims.cycle.id + add.id,
                                         median, na.rm = TRUE))
  if (sum(slices.constant) > 0)
    cat(paste(Sys.time(), ' : ', sum(slices.constant),' constant slices were found.',
            ' SSA on these for these is ommited!\n', sep=''))
                                         
  ## get slices to process                                       
  slices.process             <- !slices.constant & !slices.ocean & 
                                !slices.too.gappy & !slices.without.gaps &!slices.empty
  
  if (algorithm == 'Gapfill'){
    ##extract only a ratio of the slices to calculate for variance testing
    if (ratio.test.t != 1) {
      slices.excluded            <- logical(slices.n) 
      slices.test.n      <- ceiling(sum(slices.process) * ratio.test.t)
      ind.slices.process <- sample(which(slices.process), slices.test.n)
      slices.excluded[setdiff(which(slices.process), ind.slices.process)] <- TRUE
      slices.process[-ind.slices.process] <- FALSE
    }
    
    ##add slices usually not filled in case no validation data is available later
    if (sum(slices.process) > 0 && mean(amnt.na[slices.process]) > 0.8) {
      ind.added   <- order(amnt.na, rnorm(length(amnt.na)))[sum(slices.process)]
      slices.process[ind.added]      <- TRUE
      slices.without.gaps[ind.added] <- FALSE
    }
  }
  
  #return stuff
  if (sum(slices.process) == 0) {
    printStatus(paste('No series/slices available for filling. Most probably only',
            ' totally gappy and totally gap-free slices/series exist.', sep=''))
  }
  iters.n <- sum(slices.process)
  return(list(iters.n = iters.n, slices.process = slices.process, 
              values.constant = values.constant, slices.constant = slices.constant, 
              slices.without.gaps = slices.without.gaps, 
              slices.excluded = slices.excluded,
              slices.too.gappy = slices.too.gappy,
              tresh.min.length =  floor((1 - treshhld.gappy) * prod( dims.process.length))))
}
