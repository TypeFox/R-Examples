## written by Rick Dale (rdale@ucmerced.edu)

## build dichotomous simulated time-series data for testing
## Parameters to vary the mixing distribution of the simulated data
## values between 0,1.

#  BL1 = rand*.5; % base event rate first conditional, confederate (driver)
#  BL2 = .01; % base event rate first conditional, participant
#  BLR1 = .33; % probability of repetition
#  BLR2 = .33; % same, participant
#  BL2C1 = .25; % probability of a match, third conditional

simts <- function(BL1,BL2,BLR1,BLR2,BL2C1,tsL) {
	b1 = c(0)
	b2 = c(0)
	for (i in 1:(tsL-1)) {
	    if (runif(1) < BL1) {
	        b1 = c(b1,1)
		}
	    else if (runif(1) < BLR1 & b1[length(b1)] == 1) {
	        b1 = c(b1,1)
	    }
	    else {
	        b1 = c(b1,0)
	    }
	    
	    if (runif(1) < BL2C1 & b1[length(b1)-1] == 1) {
	    	b2 = c(b2,1)
	    }	   
	    else if (runif(1) < BL2) {
	        b2 = c(b2,1)
		}
	    else if (runif(1) < BLR2 & b2[length(b2)]==1) {
	        b2 = c(b2,1)
	    }
	    else {
	        b2 = c(b2,0)
	    }	   
	}
	return(rbind(b1,b2))	    
}


