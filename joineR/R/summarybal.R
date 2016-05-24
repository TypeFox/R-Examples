summarybal <- function(object, Y.col, times, use = "all.obs",   
    na.rm, ...)
{
    if (is.numeric(Y.col)) {
        resp <- object[, Y.col]
    }
    else {
        resp <- object[[Y.col]]
    }
    times <- unique(times)
    m.t <- apply(resp, 2, function(x){mean(x, na.rm = TRUE)})
    times <- times[order(times)]
    mean.vect <- data.frame(x = times, y = m.t)
    vv <- apply(resp, 2, function(x){var(x, na.rm = na.rm, ...)})        
    dt <- resp
    lt <- dim(resp)[2]
    cr <- matrix(ncol = lt,nrow = lt)
    diag(cr) <- 1
    
    for (i in 1 : (lt - 1)){
    	 for (j in (i + 1) : lt) {
    	      id <- !is.na( dt[, i]) & !is.na(dt[, j] )
    	      cr[i, j] <- cor(dt[id, i], dt[id, j],use = use, ...)
         }
    }
    cr[lower.tri(cr)] <- t(cr)[lower.tri(t(cr))]	    	    	
    re <- list(mean.vect = mean.vect, variance = vv, cor.mtx = cr)
    return(re)
}
