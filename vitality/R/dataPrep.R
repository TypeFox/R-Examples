#' Function for data preparation
#'
#' Function to deal with NAs, right truncated data and datatype cumulative survival or
#'  incremental motality.
#'  
#' @param time A vector of observation dates
#' @param sdata A vector of survival data of the same length as \code{time}
#' @param datatype either \code{"CUM"} for cumulative or \code{"INC"} for incremental
#' @param rc.data Boolean. Is data right-censored?
#' @param returnMatrix Boolean. False returns a data frame, true returns a matrix.
#' (as in the original), if "matrix" returns a matrix instead, with the "rc.data" column
#' being 0 for FALSE, 1 for TRUE, or 2 for TF
#' @export
#' @return Returns a data.frame or matrix with columns time, sfract, x1, x2, 
#'         Ni (incremental survival fraction), rc.data.
dataPrep <- function(time, sdata, datatype, rc.data, returnMatrix = FALSE) {
		#check for and remove NAs from data
	if (any(is.na(time))) {
        naT = is.na(time)
		time = time[!naT]
		sdata = sdata[!naT]
		warning(message = "WARNING:  NAs found in data and removed.")
	}
	if (any(is.na(sdata))) {
        naT = is.na(sdata)
		time = time[!naT]
		sdata = sdata[!naT]
		warning(message = "WARNING:  NAs found in data and removed.")
	}
	
	if(length(time) < 5) {stop(message = "ERROR:  not enough data.")}
	
	if(!all(0 <= sdata) | !all(sdata <= 1)) {
		stop("ERROR:  survival fraction data outside range.")
	}
    # end data checking	

    maxx2 <-max(time)  #for right-censored data and for plotting 
    # === check data type (CUMulative or INCremental.  If CUM, create INC ===
    if (datatype == "CUM") {
    	#====== survival assumed 1 at time 0 ===
    	if (time[1] > 0) {
    		time <- c(0, time)
    		sdata <- c(1, sdata)
    	}
    	else {
    		if (sdata[1] < 1) {
    			sdata <- sdata/sdata[1]
    			warning(message = "Initial survival < 1.  Data scaled so that initial survival = 1.")			
    		}
    	}
    
    	#------------------------
    	sfract <-sdata
    	len <- length(time)
    	#  ...set up data for MLE fitting of incremental survivorship...	
    	
    	# --------- right censored data?
    	if (rc.data != TRUE) {
    		if (rc.data == FALSE) {
    			#check if final sdata indicates full mort
    			if (sfract[len] != 0) {
    				warning("WARNING: Survival data may be right censored...")
    				rc.data<-"TF"
    			}
    			else {
    				#standard setup
    				x1 <-c(time[1:(len-1)], 0)
    				x2 <-c(time[2:len], 0)
    				sfract1 <-c(sfract[1:(len-1)], 0)
    				sfract2 <-c(sfract[2:len], 0)
    			}
    		}
    		if (rc.data == "TF") {
    			#setup: add zero sruv and short time step  ("TF" option)
    			x1 <-time
    			x2 <-c(time[2:len], 2*time[len]-time[len-1])
    			sfract1 <-sfract
    			sfract2 <-c(sfract[2:len], 0)
    		}
    		
    	}
    	else {                      #if rc.data == T
    		x1 <-time
    		x2 <-c(time[2:len], 10*maxx2)         #testing...
    		sfract1 <-sfract
    		sfract2 <-c(sfract[2:len], 0)
    	}
    	# -----------------  end of dealing with right censored data options
    	
    	Ni <- sfract1 - sfract2     #incremental survival fraction
    	
    	#  ...end conversion of cumulative survivorship to incremental mortality
    
    } else {
    	if (datatype == "INC") {
    		lent <- length(time)
    		#should be no t=0 data.  Eliminate if necessary.
    		if(time[1] == 0) {
    			time  <- time[2:lent]
    			sdata <- sdata[2:lent]
    			lent  <- length(time)
    		}
    		#check for right censored data
    		if (rc.data != T) {
    			if (rc.data == F) {
    				if(sum(sdata) < 1) {
    					rc.data <- "TF"
    					warning("WARNING: Survival data may be right censored...")
    				} else {  
    					#standard setup
    					Ni <-sdata
    					x1 <-c(0, time)[1:lent]
    					x2 <-time
    				}
    			}
    			if (rc.data == "TF") {
    				#setup: add zero sruv and short time step  ("TF" option)
    				Ni <-c(sdata, 1-sum(sdata))
    				x1 <-c(0, time)
    				x2 <-c(time, 2*time[lent]-time[lent-1])  #final time interval assumed same as prev.
    			}
    		} else {       #if rc.data == T
    			Ni <-c(sdata, 1-sum(sdata))
    			x1 <-c(0, time)
    			x2 <-c(time, 10*maxx2)    #final time interval large
    		}
    		
    		# Build cumulative data	
    	  	#sfract <-1.0      # initial survival fraction assumed 1
    	    sfract <- NULL
    	  	for (i in 1:length(Ni)) {
    			sfract <- c(sfract, 1-sum(Ni[1:i]))
    	  	}
    	    
            time <- c(0, time)
     	    len <- length(time)   #check 
    	    print(c(length(sfract), length(time), length(Ni), length(x1), length(x2), lent))
    	}
          else {
    	    stop("ERROR:  bad datatype specification")
    	}
    }
    
    if (returnMatrix) {
        returnMat <- cbind(time, sfract, Ni, x1, x2, ifelse(rc.data[1] == "TF", 2,
                                                            ifelse(rc.data[1] == "T", 1, 0)))
        dimnames(returnMat) <- list(NULL, c("time", "sfract", "Ni", "x1", "x2", "rc.data"))
        return(returnMat)    
    }
    return(data.frame(time, sfract, x1, x2, Ni, rc.data))
}
