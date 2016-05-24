#' Calculate pseudo-observations.
#' @description Computes pseudo-observations for survival analysis on right-censored data based on restricted mean survival time.
#' @usage fast_pseudo_mean(time, event, tmax)
#' @details Using a jacknife procedure and restricted mean survival time, this function calculates pseudo-observations for right-censored survival data. These pseudo-observations can be used as the response variable in a generalized estimating equations model. Missing values are not allowed in the time or event vector.
#' The function is equivalent to the pseudomean() function in the 'pseudo' package, but can handle data sets that are orders of magnitude larger.
#' @param time - Vector of follow-up times.
#' @param event - Vector of binary event statuses (0 = alive, 1 = dead).
#' @param tmax - Cut-off point for restricted mean survival time. Defaults to maximum follow-up time.
#' @examples
#' # Dummy data
#' 
#' id <- c(1, 2, 3, 4)
#' female <- c(0, 1, 1, 0)
#' time <- c(23, 45, 38, 66)
#' event <- c(1, 0, 0, 0)
#' 
#' # Compute pseudo-observations
#' 
#' pseudo = fast_pseudo_mean(time, event, 50)
#' 
#' # Create a data frame
#' 
#' test <- data.frame(id, female, pseudo)
#' 
#' # Fit a regression model
#' 
#' library(geepack)
#' 
#' summary(fit <- geese(pseudo ~ female,
#' data = test, id=id, jack = TRUE, family=gaussian,
#' corstr="independence", scale.fix=FALSE))
#' @export

fast_pseudo_mean <- function(time, event, tmax) {
  if(any(is.na(time)))
  	stop('Missing values in time vector.')
		
	if(any(time<0))
		stop('Times must be nonnegative.')
	
	if(any(is.na(event)))
		stop('Missing values in event vector.')
	
	if(any(event!=0 & event!=1))
		stop('Event must be a binary variable (alive = 0, dead = 1).')
		
	if(missing(tmax)) 
		tmax <- max(time[event==1])

	if(length(tmax)>1)
		stop('More than one value specified for tmax.')

	if(is.na(tmax))
		stop('Missing tmax value.')
  
  n <- length(time)
  df <- data.frame(id = 1:n, time, event)
  
  df$time[df$time >= tmax] <- tmax
  
  time_list <- unique(df$time)
  num_times <- length(time_list)

  summary <- do.call(
    rbind,
    lapply(
      0:tmax,
      function(x) data.frame(time = x, eligible = sum(x <= df$time), event = sum(df$event[df$time == x]))
    )
  )
  
  if(sum(summary$event[summary$time < tmax]) == 0)
  	stop('No events occured before time tmax.')
  
  n <- summary$eligible[summary['time'] == 0]
  main_rmst <- rmst_on_summary(summary, tmax)
  
  results <- do.call(
    rbind,
    mapply(
      function(x, y) {
        temp <- summary
        temp$eligible[temp$time <= x] <- temp$eligible[temp$time <= x] - 1
        
        if (y == 1) {
          temp$event[temp$time == x] <- temp$event[temp$time == x] - 1
        }
        
        this_rmst <- rmst_on_summary(temp, tmax)
        
        return(data.frame(time = x, event = y, pseudo = (n * main_rmst) - ((n - 1) * this_rmst)))
      },
      c(time_list, time_list),
      c(rep(0, num_times), rep(1, num_times)),
      SIMPLIFY = FALSE
    )
  )
  
  joined <- merge(df, results)
  
  return(joined$pseudo[order(joined$id)])
}