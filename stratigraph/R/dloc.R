dloc <- function(x, y, weights = NULL, pin.ends = FALSE,
                 method = 'genetic', end.segs = ceiling(length(x)/3),
                 pop = 100, max.gen = 200, mut = 0.01,
                 recomb = 'roulette', ext = 0.5, tol = 0.005,
                 start = 'lm', verbose = 2, plot = 1){
                 	
# dloc Draws a Line Of Correlation between two stratigraphic sections
################################################################################
# x and y give the depths of the datums in two sections. x is
#  the 'best' (i.e. reference) section
# weights is an optional vector of weights the same length as x and y
# pin.ends is TRUE to pin the loc to the first and last points; FALSE otherwise
# method, currently only 'genetic'
# end.segs is the maximum number of line segments to consider
# pop is the population size
# max.gen is the maximum number of generations to run the game
# mut is an argument multiplied by the range of possible values to
#  give the sd of the mutation applied to surviving solutions each generation
# recomb is a value between 0 and 1 giving the amount
#  of recombination among solutions that do not go extinct.
#  roulette (means fitness proportionate selection) ONLY ONE IMPLEMENTED
#  XXtournament
#  XX1 means random sampling from all living solutions
#  XX0 means no recombination; just take the best ones
# ext is the proportion of solutions that go extinct each
#  generation
# XXtol (multiplied by the sum of squared errors of the least squares linear
#  XXmodel) is the tolerance at which to end evolution of the models
# verbose and plot are flags; larger give more information
################################################################################

# to do:
#  collect _all_ lines tried
#  stopping rule (tol)
#  other weighting rules (recomb)
#  different distributions for mutation
#  keep best; wild card

################################################################################
# rationality checks

  if(length(x) == length(y)) n <- length(x)
  else stop('sections to be compared are not of same length')
  
  if(is.null(weights)) weights <- rep(1, n)
  if(length(weights) == n) weights <- as.numeric(weights)
  else stop('weights vector is the wrong length')
  
  if(sum(is.na(c(x,y))) > 0)
    stop('can not deal with NAs; remove them and try again')
  
  # if there are fewer points than end.segs; reduce the maximum
  #  number of segments so it is <= the number of points 
  if(end.segs > n){
  	warning('reducing end.segs; no point in having more segments than points')
    end.segs <- n - 1
  }
  
################################################################################
# main
  
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  
  # bestsse is the best sum of squared errors for each number of
  #  segments (segments + 1 is degrees of freedom of the model)
  bestsse <- rep(Inf, end.segs)
    
  # start with a least-squares regression line
  ls.lm <- lm(y ~ x, weights = weights)
  
  # put its sum of squared errors into the beginning of bestsse
  bestsse[1] <- sum(ls.lm$residuals^2)
  if(verbose > 0) cat('least squares regression sse:', bestsse[1], '\n')

  # scale end tolerance by starting sum of squared errors
  #tol <- bestsse[1] * tol

  # information criteria...NOT YET WORKING
  #aics <- rep(0, end.segs)
  #aics[1] <- AIC(ls.lm)
  # getAnywhere(logLik.lm)
  #aics[1] <- extractAIC(ls.lm, k = 2)[2]
  # getAnywhere(extractAIC.lm)
  
  # loc is an array with two columns giving the x and y coordinates
  #  of the vertices of the lines of correlation (which are made up of
  #  linear segments) in the living population, it is
  #  a population_size by 2 by number_of_segments array; its 'sse'
  #  attribute is a vector giving the sum of squared errors for each
  #  loc in the population
  loc <- array(rep(c(xr[1], xr[2],
                     ls.lm$fitted.values[x == xr[1]],
                     ls.lm$fitted.values[x == xr[2]]),
                   each = pop),
               dim = c(pop, 2, 2),
               dimnames = list(loc = 1:pop,
                               coord = c('x', 'y'),
                               vertex = 1:2))
  attr(loc, 'sse') <- rep(bestsse[1], pop)
  #note that at this stage all individuals in the population of lines
  #  of correleation are still the same

  # bestloc will be the best line of correlation for each number of segments
  bestloc <- vector(mode = 'list', length = end.segs)
  bestloc[[1]] <- t(loc[1,,])

  # locs will be a list of arrays giving the final population of locs for
  #  each number of segments
  locs <- vector(mode = 'list', length = end.segs)
  locs[[1]] <- loc

  if(plot){ # plot the initial linear least squares model
  	plot(x, y, main = 'initial model')
    mtext(paste('linear sse =', bestsse[1]), line = -1)
  	lines(t(bestloc[[1]]))
    points(x, ls.lm$fit, pch = 20, cex = 0.5)
    segments(x, y, x, ls.lm$fit)
    points(t(bestloc[[1]]), pch = '+')
  }

  for(segs in 2:end.segs){ # loop through different numbers of segments
  	if(verbose){
  	  cat('----------')
  	  cat(segs, 'segments')
  	  cat('----------\n')
  	}

    # initialize and fill loc array for this number of segments
    loc <- array(0, dim = c(pop, 2, segs + 1),
                 dimnames = list(loc = 1:pop,
                                 coord = c('x', 'y'),
                                 vertex = 1:(segs+1)))
    sses <- rep(0, pop)
    
    if(start == 'lm'){
      xseeds <- seq(from = xr[1], to = xr[2], length.out = segs + 1)
      yseeds <- predict(ls.lm, newdata = data.frame(x = xseeds))
      loc[,'x',] <- rep(xseeds, each = pop)
      loc[,'y',] <- rep(yseeds, each = pop)
      # add a mutation to the least squares line to initialize a population
      #  of reasonable solutions
      loc[,'x',] <- loc[,'x',] + rnorm(length(loc[,'x',]),
                           mean = 0, sd = mut * (xr[2] - xr[1]))
      loc[,'y',] <- loc[,'y',] + rnorm(length(loc[,'y',]),
                           mean = 0, sd = mut * (yr[2] - yr[1]))
    }else if(start == 'uniform'){
      for(i in 1:pop){
        xseeds <- c(xr[1], sort(runif(segs - 1, min = xr[1],
                                      max = xr[2])), xr[2])
        yseeds <- sort(runif(segs + 1, min = yr[1], max = yr[2]))    
        loc[i,'x',] <- xseeds
        loc[i,'y',] <- yseeds
      }
    }
    
    # get the sum of squared errors for all individuals in the new population
    for(ind in 1:pop){
      sse <- getsse(x, y, loc[ind,1,], loc[ind,2,], weights)
      sses[ind] <- sse$sse
    }
    attr(loc, 'sse') <- sses

    locs[[segs]] <- loc

    if(plot > 1){ # plot the starting population of solutions
      plot(x,y, main = 'starting population of locs', type = 'n')
      for(i in 1:pop) lines(t(loc[i,,]), col = 'grey')
      for(i in 1:pop) points(t(loc[i,,]), col = 'grey')
      lines(t(bestloc[[1]]))
      points(x,y)
      system('sleep 3')
    }
      
    # countdown is a counter that terminates the evolution when
    #  there has been less change than tol for a certain number of
    #  generations; initialized here
    #countdown <- 0
        
    for(gen in 1:max.gen){ # loop through generations
   	  if(verbose) cat('gen', gen, '\t')
     
      ### DO EVOLUTIONARY STUFF ###
            
      # select some locs from population (using ext)            
      fitnesses <- (max(sses) - sses) / (max(sses) - min(sses))
      fit.rank <- order(fitnesses)
      nlost <- floor(ext * length(fitnesses))

      #roulette wheel selection
      surv.index <- sort(sample(1:pop, size = pop - nlost,
                                replace = FALSE, prob = fitnesses))
      ext.index <- rep(TRUE, pop)
      ext.index[surv.index] <- FALSE
      ext.index <- (1:pop)[ext.index]
      
      # drive some extinct
      offspring <- loc
      offspring[ext.index,,] <- NA
      attr(offspring, 'sse')[ext.index] <- NA
       
      # mutate each chosen loc
      xmuts <- rnorm(pop - nlost, mean = 0, sd = mut * (max(x) - min(x)))
      ymuts <- rnorm(pop - nlost, mean = 0, sd = mut * (max(y) - min(y)))
      muts <- cbind(xmuts, ymuts)
      for(j in 1:(segs + 1)){
        offspring[surv.index,,j] <- offspring[surv.index,,j] + muts
      }
        
      # recombine some points between parents (random assortative mating)
      for(i in ext.index){
        for(j in 1:(segs + 1)){
          offspring[i,1,j] <- sample(na.omit(offspring[,1,j]), size = 1)
          offspring[i,2,j] <- sample(na.omit(offspring[,2,j]), size = 1)
        }
      }

      # constrain x-coordinates of endpoints
      offspring[,'x',1] <- xr[1] # first vertex at min(x)
      offspring[,'x',segs + 1] <- xr[2] # last vertex at max(x)
      
      # constrain y-coordinates of endpoints iff pin.ends == TRUE
      if(pin.ends == TRUE || is.character(pin.ends)){
      	if(pin.ends == TRUE || pin.ends == 'both' || pin.ends == 'bottom')
          offspring[,'y',1] <- y[x == xr[1]] # first vertex at y of smallest x
      	if(pin.ends == TRUE || pin.ends == 'both' || pin.ends == 'top')
          offspring[,'y',segs + 1] <- y[x==xr[2]] # last vertex at y of largest x
      }

      # constrain each loc to be monotonic non-decreasing in x and y
      for(i in 1:pop){
        offspring[i,'y',] <- sort(offspring[i,'y',])
        offspring[i,'x',] <- sort(offspring[i,'x',])
      }
      # note: occasionally this will swap points; cf. gene conversions
            
      # [peak climb best solution?]

      ### END EVOLUTIONARY STUFF ###

      if(plot > 1){
      	plot(x,y, type = 'n', main = paste(segs, 'segments; generation', gen))
        for(i in surv.index) lines(t(loc[i,,]), col = 'green')
        for(i in ext.index) lines(t(loc[i,,]), col = 'red')
        #if(!is.matrix(loc[fitnesses == 1, , ])) browser()
        #lines(t(loc[fitnesses == 1,,]), col = 'grey')
        # this line seems to die with: arg. to t() is not a matrix
        points(x,y, pch = 20)
        mtext('survivors green; extinct solutions red; best parent grey')
        system('sleep 0.5')
        
        plot(x,y, type = 'n', main = paste(segs, 'segments; generation', gen))
      	for(i in surv.index) lines(t(offspring[i,,]), col = 'blue')
      	for(i in ext.index) lines(t(offspring[i,,]))
      	points(x,y, pch = 20)
        mtext('mutations blue; offspring black')
        system('sleep 0.5')
      }

      # update loc and recalculate sses
      loc <- offspring
      for(i in 1:pop){
        sse <- getsse(x, y, loc[i,1,], loc[i,2,], weights)
        sses[i] <- sum(sse$res^2, na.rm = TRUE)
      }
      attr(loc, 'sse') <- sses
  	  
  	  thisbestsse <- min(sses)
      bestone <- (1:pop)[sses == min(sses)][1]
  	  thisbestloc <- t(loc[bestone,,])

      # if the decrease in sse is less than tol; count down...
      #if(abs(bestsse[segs] - thisbestsse) < tol){
      #  countdown <- countdown + 1
      #}else{
      #	 #otherwise reset countdown
      #  countdown <- 1
      #}

   	  if(thisbestsse < bestsse[segs]){
   	  	bestsse[segs] <- thisbestsse
   	  	bestloc[[segs]] <- thisbestloc
   	  }
   	     	  
      if(verbose){
      	cat('thisbestsse:', thisbestsse, '\t')
        cat('bestsse:', bestsse[segs], '\n')
#        cat('tol', countdown, '\n')
      }
      
      # after 10 generations with no decrease smaller than tol
      #if(countdown > 10) break
      } # end looping through generations
    
    # print best sse of this number of segments
    if(verbose) cat('best sse with', segs, 'segments:', bestsse[segs], '\n')
    
    # plot best loc of this number of segments
    #if(plot){
    #  points(x, sse$fit, col = 'red')
    #  lines(sse$mod, col = 'green')
    #  segments(x, y, x, sse$fit)
    #}
    
    # add best loc to locs list
    locs[[segs]] <- loc
  } # end loop for this number of segments
  
  # produce a final summary plot
  if(plot){
    plot(x, y, main = 'best loc for each number of segments', type = 'n')
    for(i in seq_along(locs)){
      if(is.null(locs[[i]])) next
      for(j in 1:pop){
      	if(i == 1) lines(locs[[i]][j,,], col = 'grey')
      	else lines(t(locs[[i]][j,,]), col = 'grey')
      }
    }
    for(i in seq_along(locs)){
    	if(i == 1) lines(t(bestloc[[i]]), col = i)
    	else lines(bestloc[[i]], col = i)
    }
    points(x,y)
    legend('topleft', legend = seq_along(locs),
            lwd = 1, col = seq_along(locs))
  }
  return(list(locs = locs, bestloc = bestloc, bestsse = bestsse))
}

################################################################################
# sub

## getsse
# given raw data in x and y and a model (loc) described by two vectors,
#  xvert and yvert, giving the x and y coordinates of its vertices,
#  return the model, fit, residuals and sum of squared errors 
getsse <- function(x, y, xvert, yvert, weights){
  res <- rep(0, length(x))
  fit <- rep(0, length(x))
  for (i in 1:length(x)){ # for each point...
  	if(all(x[i] < xvert)){ # extrapolate first two points
  	  x2 <- xvert[2]
      x1 <- xvert[1]
      y2 <- yvert[2]
      y1 <- yvert[1]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1 
      fit[i] <- m * x[i] + b
  	}else if(all(x[i] > xvert)){ # extrapolate last two points
  	  x2 <- xvert[length(xvert)]
      x1 <- xvert[length(xvert) - 1]
      y2 <- yvert[length(xvert)]
      y1 <- yvert[length(xvert) - 1]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1 
      fit[i] <- m * x[i] + b
  	}else if(x[i] %in% xvert){ # if model matches data exactly
  	  fit[i] <- NA
  	}else{ # interpolate points inside x-range of model
      smaller <- max(which(xvert < x[i])) # vertex left of ith point
      larger <- min(which(xvert > x[i])) # vertex right of ith point
      x2 <- xvert[larger]
      x1 <- xvert[smaller]
      y2 <- yvert[larger]
      y1 <- yvert[smaller]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1 
      fit[i] <- m * x[i] + b
    }
  }
  fit[is.na(fit)] <- yvert[xvert %in% x]
  mod <- list(x = xvert, y = yvert)
  res <- y - fit
  sse <- sum((weights * res)^2)
  return(list(mod = mod, fit = fit, res = res, sse = sse))
}

## cpredict
# given a model (loc) described by two vectors, xvert and yvert.
#  and new data in x, return the predicted y values 
cpredict <- function(xvert, yvert, newx){
  if(length(newx) == 0) return(numeric(0))
  newy <- rep(0, length(newx))
  for (i in 1:length(newx)){ # for new x value...
  	if(all(newx[i] < xvert)){ # extrapolate first two points
  	  x2 <- xvert[2]
      x1 <- xvert[1]
      y2 <- yvert[2]
      y1 <- yvert[1]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1 
      newy[i] <- m * newx[i] + b
  	}else if(all(newx[i] > xvert)){ # extrapolate last two points
  	  x2 <- xvert[length(xvert)]
      x1 <- xvert[length(xvert) - 1]
      y2 <- yvert[length(xvert)]
      y1 <- yvert[length(xvert) - 1]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1 
      newy[i] <- m * newx[i] + b
  	}else if(newx[i] %in% xvert){ # if model matches data exactly
  	  newy[i] <- NA
  	}else{ # interpolate points inside x-range of model
      smaller <- max(which(xvert < newx[i])) # vertex left of ith point
      larger <- min(which(xvert > newx[i])) # vertex right of ith point
      x2 <- xvert[larger]
      x1 <- xvert[smaller]
      y2 <- yvert[larger]
      y1 <- yvert[smaller]
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1
      newy[i] <- m * newx[i] + b
    }
  }
  newy[is.na(newy)] <- yvert[xvert %in% newx]
  return(newy)
}