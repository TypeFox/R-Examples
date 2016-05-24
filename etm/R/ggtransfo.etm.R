### Function to transform a etm object into something usable for ggplot ###
### Arthur Allignol <arthur.allignol@uni-ulm.de                         ###

ggtransfo <- function(x, ...) {
    UseMethod("ggtransfo")
}

ggtransfo.etm <- function(x, tr.choice, ...) {

    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")

    sx <- summary(x, ...)
    
    if (missing(tr.choice)) tr.choice <- names(sx)
    sx <- sx[tr.choice]
    
    sx_display <- do.call(rbind, lapply(seq_along(tr.choice), function(i) {
        tmp <- sx[[i]]
        tmp$trans <- tr.choice[i]
        
        tmp$timemax <- c(tmp$time[-1], max(tmp$time) + 1)
        tmp
    }))
    
    sx_display
}


### test
## aa <- ggtransfo(tr.prob, "0 1")

## p <- ggplot(aa) + 
##     geom_step(aes(x = time, y = P)) +
##     geom_rect(aes(xmin = time, xmax = timemax, ymin = lower, ymax = upper),
##               alpha = 0.3) 
