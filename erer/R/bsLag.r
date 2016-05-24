bsLag <- function(h, lag, prefix = "", var.name,
    suffix = ".t_", include.orig = TRUE, by.lag = FALSE, ...)
{
    if(!is.ts(h) ) {stop("Please provide time series data.\n")}
    mh <- as.matrix(h); n1 <- dim(mh)[1]; n2 <- dim(mh)[2]
    if(lag >= n1) {stop("Lag should be less than data dimension.\n")}
    
    if (missing(var.name)) {
        if (n2==1) {var.name <- deparse(substitute(h))
        } else {var.name <- colnames(h)}
    } else {
        if(length(var.name)!=n2 ) {
        stop("Length of 'var.name' should equal to variable numbers")}
    }

    all <- data.frame(matrix(nrow = n1 - lag, ncol = n2 * (lag + 1)))  
    for (i in 1:n2) {       
      out <- data.frame(embed(mh[,i], dimension=lag+1) ) 
      ww <- paste(prefix, var.name[i], suffix, "0", sep="")
      if (lag >= 1) {
          for(j in 1:lag) {
              ww2 <- paste(prefix, var.name[i], suffix, j, sep="")
              ww <- c(ww, ww2)
          }
      }        
      sel <- if (by.lag) {
        seq(from = i, to = n2 * lag + i, by = n2 )  
      } else {
        seq(from = (lag + 1) * (i - 1) + 1, to = (lag + 1) * i, by = 1)
      }
      all[, sel] <- out                   
      colnames(all)[sel] <- ww
    }

    del <- grep(paste(suffix, '0', sep=''), colnames(all))
    out <- if(include.orig) {all} else {all[, -del ]}   
    fia <- ts(out, start=start(h) + c(0, lag), end=end(h), 
      frequency=tsp(h)[3])
    return(fia)
} 