# ---Author: Jake Zwart, 2013-09-30 --- 
# cleaned from previous script - original author: Arianto Santoso 

#'@export
k.cole <- function(ts.data){
	if(!has.vars(ts.data, 'wnd')){
		stop('k.cole requires a "wnd" column in the supplied data')
	}
	wnd <- get.vars(ts.data, 'wnd')
	k600 <- k.cole.base(wnd[,2])

	return(data.frame(datetime=ts.data$datetime, k600=k600))
}


# wnd: wind value in m/s

# OUTPUT: returns the gas exchange velocity for k600 in units of m/day
#'@export
k.cole.base <- function(wnd){
  U10 <- wnd  #This function uses just the wind speed it is supplied. 
  k600 <- 2.07 + (0.215 * (U10^1.7)) # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}

