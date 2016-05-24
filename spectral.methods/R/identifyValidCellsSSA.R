.identifyValidCellsSSA = function(
##title<< helper function for Ncdf SSA algorithms 
    dims.cycle.id, dims.process.id, datacube, ratio.const, tresh.const , print.status,
    slices.n, algorithm,  g = c() ,process.cells = c('gappy','all')[1], dims.cycle = c(), 
    args.call.SSA = list(), tresh.fill.dc = 0, ratio.test.t =1, first.guess = 'mean', 
    ocean.mask = c(), dims.process = c(), dims.process.length = 0, 
    slices.without.gaps = rep(FALSE, slices.n), slices.too.gappy = rep(FALSE, slices.n),
    slices.constant = rep(FALSE, slices.n), slices.process = rep(TRUE, slices.n),
    slices.ocean = rep(FALSE, slices.n), values.constant = integer(length = slices.n),
    slices.excluded = rep(FALSE, slices.n), MSSA = FALSE)   
##details<< helper function for gapfillNCdf and decomposeNcdf that identifies the grid cells to process. 
##seealso<<
##\code{\link{gapfillNCdf}}, \code{\link{gapfillNCdf}}    
{  
  ## TODO
  ## -possibility to identify gap less MSSA blocks
  ## -include possibility to infer slices.continuous max borders

  if (print.status)
    cat(paste(Sys.time(), ' : Identifying valid cells ...\n', sep=''))
  
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
          n_biggaps   <- max(c(floor(dtpts * gaps[1] /size.bg), 1))  
          if (gaps[1] > 0) {
            n_smallgaps <- n_biggaps * size.bg * gaps[1] / gaps[2]
          } else {
            n_smallgaps <- dtpts * gaps[2]
          }
          treshhld.gappy <- 1 - ((n_smallgaps + n_biggaps * size.bg) / dtpts ) - tresh.fill.dc
        } else {
          treshhld.gappy <- 1 - (tresh.fill.dc)
        }
        slices.too.gappy <- as.vector(amnt.na >= treshhld.gappy)
        
        ## modify this in case first guess is supplied
        if ((length(first.guess) > 1) & tresh.fill.dc == 0) {
          amnt.na.first.guess    <- apply(first.guess, MARGIN = dims.cycle.id + 1,
              function(x)sum(is.na(x)) / dtpts  )
          if ((length(ocean.mask) > 0 ) & (sum(!is.na(match(c('longitude','latitude'), dims.process))) == 2))
            amnt.na.first.guess  <- 1- apply(first.guess, MARGIN = dims.cycle.id + 1 ,
                function(x) sum(!is.na(x[as.vector(!ocean.mask)])) / sum(!ocean.mask)   )
          slices.too.gappy[amnt.na.first.guess < 0.9 & amnt.na < 0.75] <- FALSE
        }      
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
              slices.too.gappy = slices.too.gappy))
}
