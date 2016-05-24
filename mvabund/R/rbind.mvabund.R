###############################################################################
## rbind.mvabund: combine two or more mvabund object columnwise to another 	  #
## mvabund object											                                        #
###############################################################################
rbind.mvabund <- function(..., deparse.level = 1) {
objects <- list(...)

objects <- lapply(objects, function(x) {
       unabund(x)
    })

mrbind <-  do.call( "rbind", c(objects ,deparse.level=deparse.level))

return(as.mvabund(mrbind))

}


