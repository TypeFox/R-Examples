
################################################
# Computation of a progress bar
# INPUT:
# ops ... number of operations (loop index)
# prblen ... length of progress bar
#---------------------------------
BIFIE.progressbar <- function( ops , prblen ){
    prblen -> prb
    vec <- seq( 1 , ops  )
    vec[ ops ] <- ops - .1
    NR <- ops / prb
    m1 <- vec %% NR
    pr1 <- 1 * ( diff(m1) < 0  )
    pr1 <- c( 1 , pr1 )   
	# returns a vector of zeroes and one indicating
	# iteration of a move in th progress bar
    return(pr1)
        }
####################################################