
###############################################################################
## cbind.mvabund: combine two or more mvabund object columnwise to another 	  #
## mvabund object											                                        #
###############################################################################
cbind.mvabund <- function(..., deparse.level = 1) {
objects <- list(...)

objects <- lapply(objects, function(x) {
       unabund(x)
    })

mcbind <-  do.call( "cbind", c(objects ,deparse.level=deparse.level))

return(as.mvabund(mcbind))

}

