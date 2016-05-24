plot.catIrt <-
function(x, which = "all", ids = "none", conf.lev = .95, legend = TRUE, ask = TRUE, ... ){
 
# First, making sure that if "ids" or "which" is all, we should include all plots.
  if( any(which == "all") )
    which <- 1:4
  	
  if( any(ids == "all") )
    ids <- seq_along(x$cat_theta)
  	
# Next, checking the "which" and "ids" argument.
  if( !all(which == 1 | which == 2 | which == 3 | which == 4) & any(which != "none") )
  	stop( "which must be 'all', 'none', or an integer vector containing 1 -- 4" )
  	
  if( any(ids != "none" & ids != "all") )
    if( missing(ids) | ( any( ids < 1 | ids > length(x$cat_theta) ) ) )
      stop( "ids must be 'all', 'none', or an integer vector containing 1 -- length(theta)" )
  
# Making sure the conf.lev is within the bounds (only if there are ids to plot)).
  if( all( ids != "none" ) ){
    if( conf.lev < 0 | conf.lev > 1 ){
  	  stop("conf.lev must be a positive number between 0 and 1")
  	}
  } # END if STATEMENT


##########################
# Total Bank Information #
##########################

  info.t <- vector(mode = "list", length = 3)
  info.e <- vector(mode = "list", length = 2)
  info.o <- vector(mode = "list", length = 2)
  
# Figuring out the realistic minimum and maximum of our data.

# First: find the min and max of the theta.
  min_theta <- min( unique(min(x$tot_theta)), unique(min(x$full_params[ , -1])) )
  max_theta <- max( unique(max(x$tot_theta)), unique(max(x$full_params[ , -1])) )
  
# Second: build a sequence from the lowest to the highest.
  info.t$eval <- seq(min_theta, max_theta, length.out = 1000)
  
# Third: pull out parameters and indicate class of those parameters.
  full_params        <- x$full_params[ , -c(1, ncol(x$full_params))]
  class(full_params) <- c(x$mod$mod, "matrix")

# Fourth: calculate info on those parameters across all people.
  y            <- FI(params = full_params, theta  = info.t$eval, type = "expected")
  info.t$info  <- y$test
  info.t$sem   <- y$sem
  
  for( i in seq_along(x$cat_theta) ){
  	
# Pulling out parameters for person 'i'.
    cat_params        <- x$cat_indiv[[i]]$cat_params[ , -c(1, ncol(x$full_params))]
    class(cat_params) <- c(x$mod$mod, "matrix")
  	
# Calculating observed and expected info on that person.
    z1 <- FI(params = cat_params,
             theta  = x$cat_theta[i],
             type   = "expected" )
    z2 <- FI(params = cat_params,
  	         theta  = x$cat_theta[i],
             type   = "observed",
             resp   = x$cat_indiv[[i]]$cat_resp )
  	
# Storing observed and expected information.
  	info.e$info[i] <- z1$test
  	info.o$info[i] <- z2$test
  	info.e$sem[i]  <- z1$sem
  	info.o$sem[i]  <- z2$sem
  }
  

##################
# Begin Plotting #
##################
  
  if( ask ){
    eask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(eask))
  }
  
########################
# Plotting Information #
########################
	
# Plotting the Bank Information Function.
  if( any(which == 1) ){
    plot(x = info.t$eval, y = info.t$info,
         type = "l", lwd = 2,
         main = "Bank Information Function",
         xlab = expression(theta), ylab = "Expected Information", axes = FALSE)
    axis(1, col = "grey")
    axis(2, col = "grey", las = 1)
  } # END if STATEMENT

# Plotting the Bank SEM Function.
  if( any(which == 2) ){
    plot(x = info.t$eval, y = info.t$sem,
         type = "l", lwd = 2,
         main = "Bank SEM Function",
         xlab = expression(theta), ylab = "Expected SEM", axes = FALSE)
    axis(1, col = "grey")
    axis(2, col = "grey", las = 1)
  } # END if STATEMENT

# Plotting the Empirical and Observed Information Functions.
  if( any(which == 3) ){
    plot(x = sort(x$cat_theta), y = info.e$info[order(x$cat_theta)],
         type = "l", lty = 1, lwd = 2,
         main = "Obs and Expected CAT Fisher Information",
         xlab = expression(theta), ylab = "Information", axes = FALSE)
    axis(1, col = "grey")
    axis(2, col = "grey", las = 1)
    lines(x = sort(x$cat_theta), y = info.o$info[order(x$cat_theta)],
          type = "l", lty = 2, lwd = 2, col = 2 )
          
    if( legend ){     
      legend(x = "bottomright",
             y = c("Expected", "Observed"),
             lty = 1:2, lwd = 2, col = 1:2)
    } # END if STATEMENT
    
  } # END if STATEMENT
  
# Plotting the Empirical and Observed SEM Functions.
  if( any(which == 4) ){
    plot(x = sort(x$cat_theta), y = info.e$sem[order(x$cat_theta)],
         type = "l", lty = 1, lwd = 2,
         main = "Obs and Expected CAT Fisher-Based SEM",
         xlab = expression(theta), ylab = "SEM", axes = FALSE)
    axis(1, col = "grey")
    axis(2, col = "grey", las = 1)
    lines(x = sort(x$cat_theta), y = info.o$sem[order(x$cat_theta)],
          type = "l", lty = 3, lwd = 2, col = 3 )
    
    if( legend ){    
      legend(x = "topright",
             y = c("Expected", "Observed"),
             lty = c(1, 3), lwd = 2, col = c(1, 3))
    } # END if STATEMENT
     
  } # END if STATEMENT

  
###################################
# Plotting Individual Progression #
###################################

# Plotting CAT Steps for Desired People (if it is requested).
  if( !is.null(ids) & any(ids != "none") ){

# For all people do the following:
# - pull out the id number,
# - pull out the vector of theta estimates for that person in the cat,
# - plot the theta estimates along with 1.96 times their standard errors,
# - write nice things, and include a line for the Bank Ability estimate.
  	for( i in seq_along(ids) ){
  	  z.crit                <- qnorm( p = (1 - conf.lev)/2 )
      id.i                  <- ids[i]
      theta.i               <- x$cat_indiv[[id.i]]$cat_theta
      sem.i                 <- x$cat_indiv[[id.i]]$cat_sem
      sem.i[ is.na(sem.i) ] <- 10
      
# Build a list to make calling plots more customizable
      p.list <- list(xlab = "CAT Step", ylab = expression(theta),
                     main = paste( "Ability Estimates for Simulee ", id.i, sep = ""),
                     type = "b", lty = 2, lwd = 2, pch = 20, axes = FALSE)
      p.list <- modifyList(p.list, list(x = seq_along(theta.i), y = theta.i, ...) )
        
# Plot several things, including one customizable and the others not-so-much
      do.call(plot, p.list)
      
      if( !p.list$axes ){
        axis(1, col = "grey",
             at = seq_along(theta.i), labels = seq_along(theta.i) - 1)
        axis(2, col = "grey", las = 1)
      } # END if STATEMENT
      
      mtext(text =  paste( "CAT Theta Estimate: ",  round(x$cat_theta[id.i], 3), "; ",
                           "Bank Theta Estimate: ", round(x$tot_theta[id.i], 3), sep = "" ),
            side = 3, line = .5)
      segments(x0 = seq_along(theta.i), x1 = seq_along(theta.i),
               y0 = theta.i - z.crit * sem.i,
               y1 = theta.i + z.crit * sem.i,
               col = "red", lty = 3 )
      abline(h = x$tot_theta[id.i],
             lty = 4, lwd = 2, col = "burlywood4")
      
# Printing the bounds, adding true theta lines, and making a legend:
      y.leg <- c.leg <- l.leg <- w.leg <- NULL

#####
# 1 # (IF TRUE THETA)
#####
      if( length(x$true_theta) >= 1 & all( !is.na(x$true_theta) ) ){      	
        abline(h = x$true_theta[id.i],
               col = "blue", lty = 5, lwd = 2)
        y.leg <- c(y.leg, paste("True Theta = ", round(x$true_theta[id.i], 3), sep = ""))
        c.leg <- c(c.leg, "blue")
        l.leg <- c(l.leg, 5)
        w.leg <- c(w.leg, 2)      	
      } # END if STATEMENT
      
#####
# 2 # (BANK ABILITY ESTIMATE)
#####
      y.leg <- c(y.leg, "Bank Theta Estimate")
      c.leg <- c(c.leg, "burlywood4")
      l.leg <- c(l.leg, 4)
      w.leg <- c(w.leg, 2)
      
#####
# 3 # (CLASSIFICATION BOUNDS)
#####
      if( any(x$mod$catTerm$term == "class") ){
        abline(h = x$mod$catTerm$c.term$bounds[id.i, ],
               col = "green")
        y.leg <- c(y.leg, "Classification Bound(s)")
        c.leg <- c(c.leg, "green")
        l.leg <- c(l.leg, 1)
        w.leg <- c(w.leg, 1)
      } # END if STATEMENT
      
#####
# 4 # (CONFIDENCE INTERVAL)
#####
      if( conf.lev > 0 ){
        y.leg <- c(y.leg, "Confidence Interval")
        c.leg <- c(c.leg, "red")
        l.leg <- c(l.leg, 2)
        w.leg <- c(w.leg, 1)
      } # END if STATEMENT

#####
# 5 # (PRINTING THE LEGEND)
#####
      if( legend ){
        legend( x = "topright",
                y = y.leg, col = c.leg, lty = l.leg, lwd = w.leg )
      } # END if STATEMENT
      
    } # END for i LOOP
  } # END if STATEMENT
 

####################
# Returning Values #
####################

# And returning values so people can plot for themselves:  
  ans <- list(bank_info = list( eval = sort(info.t$eval),
                                info = info.t$info[order(info.t$eval)],
                                sem  = info.t$sem[order(info.t$eval)] ),                   
              exp_info  = list( eval = sort(x$cat_theta),
                                info = info.e$info[order(x$cat_theta)],
                                sem  = info.e$sem[order(x$cat_theta)] ),
              obs_info  = list( eval = sort(x$cat_theta),
                                info = info.o$info[order(x$cat_theta)],
                                sem  = info.o$sem[order(x$cat_theta)] ) )
                                
  invisible(ans)
  
} # END plot.catIrt FUNCTION