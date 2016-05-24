## Making dBMvar a generic funtion if (!isGeneric('dBMvar')) {
## setGeneric('dBMvar', function(BMvars, BMvar, n.locs, break.list)
## standardGeneric('dBMvar')) } if (!isGeneric('brownian.motion.variance.dyn'))
            brownian.motion.variance <- function(time.lag, location.error, x, y) {
              # Creating NULL vectors to store data
              n.locs <- unique(c(length(time.lag), length(location.error), length(x), length(y)))
              if (length(n.locs) != 1) 
                stop("Not an equal number of locations in call to brownian.motion.variance")
              T.jump <- alpha <- ztz <- loc.error.1 <- loc.error.2 <- NULL
              i <- 2
              while (i < n.locs) {
                # calculate serveral variables needed in the variance function
                t <- time.lag[i - 1] + time.lag[i]
                T.jump <- c(T.jump, t)
                a <- time.lag[i - 1]/t
                alpha <- c(alpha, a)
                u <- c(x[i - 1], y[i - 1]) + a * (c(x[i + 1], y[i + 1]) - c(x[i - 1], 
                                                                            y[i - 1]))
                ztz <- c(ztz, (c(x[i], y[i]) - u) %*% (c(x[i], y[i]) - u))
                loc.error.1 <- c(loc.error.1, location.error[i - 1])
                loc.error.2 <- c(loc.error.2, location.error[i + 1])
                i <- i + 2
              }
              # Likelihood function for Brownian Motion variance estimation
              likelihood <- function(var, T.jump, alpha, loc.error.1, loc.error.2, ztz) {
                v <- T.jump * alpha * (1 - alpha) * var + ((1 - alpha)^2) * (loc.error.1^2) + 
                  (alpha^2) * (loc.error.2^2)
                l <- log((1/(2 * pi * v)) * exp(-ztz/(2 * v)))
                if (any(is.na(l))) 
                  stop("An internal error occoured. Contact maintainer.")
                return(-sum((l)))
                # sum was using na.rm=T, but i want to know why probably has to do with not
                # possible optimizations, now build in logical statement in line before
              }
              BMvar <- optimize(likelihood, lower = (l <- 0), upper = (u <- 1e+15), T.jump = T.jump, 
                                alpha = alpha, loc.error.1 = loc.error.1, loc.error.2 = loc.error.2, 
                                ztz = ztz)  # implement checks if optimization worked
              if (any(BMvar$minimum %in% c(l, u))) 
                stop("Optimization failed! Consider changing mapunits")
              
              if ((length(x)%%2) != 1) 
                warning("Not an even number of locations in variance function")
              
              return(list(BMvar = BMvar$minimum, cll = -BMvar$objective))
            }
setGeneric("brownian.motion.variance.dyn", function(object, location.error, window.size, margin) {standardGeneric("brownian.motion.variance.dyn")})

setMethod(f = "brownian.motion.variance.dyn", 
          signature = c(object = ".MoveTrackSingle", location.error = "numeric", window.size = "numeric", margin = "numeric"), 
          definition = function(object, location.error, window.size, margin) {
            time.lag <- timeLag(object, units = "mins")  #units need to correspont between BBMM method and here
	  if(n.locs(object)!= length(location.error))
		  stop("The location error vector has not the same length as the move object")
	  if(n.locs(object)<window.size)
		  stop("window.size can't be larger than the number of locations in the move object")
	  if(any(is.na(location.error)))
		  stop("The location error contains NAs")
            if(isLonLat(object)) stop("You can not use longitude latitude projection for this function. To transform your coordinates use the spTransform function. \n")
            if (any((c(margin, window.size)%%2) != 1)) 
              stop("Margin and window size need to be odd")
            # function to calculate brownian.motion.variance for a piece of track
            
            breaks <- margin:(window.size - margin + 1)
	        if(is.unsorted(breaks)){
			      stop("The proposed breaks are unsorted, is the window.size smaller than 2*margin?")
	        }
            uneven.breaks <- breaks[(breaks%%2) == 1]
            breaks.found <- c()
            BMvars <- data.frame(BMvar = c(), loc = c())
            
            if (length(breaks) < 2) 
              stop("Margin to window ratio not appropriate")
            # What would be necessary to change?  try all possible window starts in track
            for (w in 1:(n.locs(object) - window.size + 1)) {
              # calculate all vectors for parts inside the window
              x.sub <- coordinates(object)[w:(w - 1 + window.size), 1]
              y.sub <- coordinates(object)[w:(w - 1 + window.size), 2]
              time.lag.sub <- time.lag[w:(w - 1 + window.size)]
              location.error.sub <- location.error[w:(w - 1 + window.size)]
              # calculate bmvar for whole window
              wholeWindow <- brownian.motion.variance(x = x.sub, y = y.sub, location.error = location.error.sub, 
                                                      time.lag = time.lag.sub)
              wholeWindow$BIC <- -2 * wholeWindow$cll + log(window.size)
              
              breakWindow <- list(BIC = Inf)
              ## is a list use inf so all other numbers are lower try all possible breaks and
              ## find the one with minimal BIC, b represents break location
              for (b in uneven.breaks) {
                before <- brownian.motion.variance(x = x.sub[1:b], y = y.sub[1:b], location.error = location.error.sub[1:b], 
                                                   time.lag = time.lag.sub[1:b])
                after <- brownian.motion.variance(x = x.sub[b:window.size], y = y.sub[b:window.size], 
                                                  location.error = location.error.sub[b:window.size], time.lag = time.lag.sub[b:window.size])
                # minimize aic check N (window.size) here TODO
                breakBIC <- (-2 * (before$cll + after$cll) + 2 * log(window.size))
                if (breakBIC < breakWindow$BIC) 
                  breakWindow <- list(w = w, b = b, var.before = before$BMvar, var.after = after$BMvar, 
                                      BIC = breakBIC)
              }
              # breakWindow should now be the best possible break
              
              if (breakWindow$BIC < wholeWindow$BIC) {
                # check if the break is better than any other, if so use those variance values
                windowBMvar <- c(rep(breakWindow$var.before, sum(breaks < breakWindow$b)), 
                                 rep(breakWindow$var.after, sum(breaks > breakWindow$b)))  # use > instead of >= because the sigma is going to be assosicated with the segments not the locations and there are one more locations than segments in the track
                breaks.found <- c(breaks.found, (w - 1 + breakWindow$b))
              } else {
                windowBMvar <- rep(wholeWindow$BMvar, length(breaks) - 1)
              }
              BMvars <- rbind(BMvars, data.frame(BMvar = windowBMvar, loc = w - 1 + margin:(window.size - 
                margin)))
            }
            
            tmp <- aggregate(BMvar ~ loc, data = BMvars, function(x) {
              c(mean = mean(x), length = length(x))
            })
            if (is.null(breaks.found)) 
              breaks.found <- numeric()
            DBMvar <- new("dBMvariance", 
                          object,
                          margin = margin, 
                          window.size = window.size, 
                          means = c(rep(NA, min(tmp$loc) - 1), tmp$BMvar[, "mean"], rep(NA, n.locs(object) - max(tmp$loc))), 
                          in.windows = c(rep(NA, min(tmp$loc) - 1), tmp$BMvar[, "length"], rep(NA, n.locs(object) - max(tmp$loc))), 
                          interest = c(rep(FALSE, min(tmp$loc) - 1), tmp$BMvar[, "length"] == max(tmp$BMvar[, "length"]), rep(FALSE, n.locs(object) - max(tmp$loc))), 
                          break.list = breaks.found)
            return(DBMvar)
          })


setMethod(f = "brownian.motion.variance.dyn", signature = c(object = ".MoveTrackSingleBurst", 
							    location.error = "numeric", window.size = "numeric", margin = "numeric"), definition = function(object, 
							    location.error, window.size, margin) {
	res <- brownian.motion.variance.dyn(as(object, ".MoveTrackSingle"), location.error = location.error, 
					    window.size = window.size, margin = margin)
	return(new("dBMvarianceBurst", as(res, "dBMvarianceTmp"), object))
	 })
