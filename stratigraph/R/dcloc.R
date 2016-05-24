dcloc <- function(x, y, pins, weights = NULL,
                 method = 'genetic', end.segs = 3,
                 pop = 100, max.gen = 200, mut = 0.01,
                 recomb = 'roulette', ext = 0.5, tol = 0.005,
                 start = 'lm', verbose = 2, plot = 1){

# dcloc Draws a Constrained Line Of Correlation between two stratigraphic sections
################################################################################
# rationality checks

  if(length(x) == length(y)) n <- length(x)
  else stop('sections to be compared are not of same length')
  
  if(is.null(weights)) weights <- rep(1, n)
  if(length(weights) == n) weights <- as.numeric(weights)
  else stop('weights vector is the wrong length')
  
  if(sum(is.na(c(x,y))) > 0)
    stop('can not deal with NAs; remove them and try again')
  
  if(is.null(pins)){
    stop('use dloc() if there are no constraints')
  }else{
    if(length(pins$x) != length(pins$y))
      warning('x and y coordinates for pins different lengths')
    pins <- as.data.frame(lapply(pins, as.numeric))
    if(any(pins$x != sort(pins$x))){
      warning('pins out of order; reordering it by x')
      preorder <- order(pins$x)
      pins$x <- pins$x[preorder]
      pins$y <- pins$y[preorder]
      if(any(pins$y != sort(pins$y)))
        stop('y coordinates of pins not monotonically increasing with x')
    }
    cat(nrow(pins), 'pins;', nrow(pins) + 1, 'independent solutions\n')
  }

################################################################################
# main
  
  dl <- vector(mode = 'list', length = nrow(pins) + 1)
  bestloc <- vector(mode = 'list', length = end.segs) 
  for(p in 1:(nrow(pins) + 1)){

    xlt <- pins[p,'x']
    xgt <- pins[p - 1,'x']
    if(length(xgt) == 0) xgt <- -Inf
    if(is.na(xlt)) xlt <- Inf
    wh <- (x < xlt) & (x > xgt)

  	if(p == 1){ # first independent solution
      thisx <- c(x[wh], pins[1,'x'])
      thisy <- c(y[wh], pins[1,'y'])
      thisweights <- c(weights[wh], 1)
      pin.ends <- 'top'
    }else if(p == nrow(pins) + 1){ #last independent solution
      thisx <- c(pins[nrow(pins),'x'], x[wh])
      thisy <- c(pins[nrow(pins),'y'], y[wh])
      thisweights <- c(1, weights[wh])
      pin.ends <- 'bottom'
    }else{ #solutions between pins
      thisx <- c(pins[p - 1,'x'], x[wh], pins[p,'x'])
      thisy <- c(pins[p - 1,'y'], y[wh], pins[p,'y'])
      thisweights <- c(1, weights[wh], 1)
      pin.ends <- TRUE
    }
    # if there are fewer points than end.segs; reduce the maximum
    #  number of segments so it is <= the number of points 
    if(end.segs > length(thisx)){
      warning('reducing end.segs; no point in having more segments than points')
      end.segs <- end.segs - 1
    }
    if(verbose > 1){
      cat('--------------------')
      cat('independent solution', p, ';', length(thisx), 'points')
      cat('--------------------\n')
    }
        
    dl[[p]] <- dloc(thisx, thisy, thisweights, pin.ends = pin.ends,
                    method = method, end.segs = end.segs,
                    pop = pop, max.gen = max.gen, mut = mut,
                    recomb = recomb, ext = ext, tol = tol,
                    start = start, verbose = 2, plot = 1)
  }
  
  for(p in 1:(nrow(pins) + 1)){
    for(segs in 2:end.segs){
      bestloc[[segs]] <- rbind(bestloc[[segs]], dl[[p]]$bestloc[[segs]])
    }
  }

  if(plot > 1){
    plot(x,y)
    points(pins, pch = 19, col = 'blue')
    for(i in seq_along(bestloc)) lines(bestloc[[i]], col = i)
    legend('topleft', legend = seq_along(bestloc),
           lwd = 1, col = seq_along(bestloc))

  }
  return(list(dl = dl, bestloc = bestloc))
}