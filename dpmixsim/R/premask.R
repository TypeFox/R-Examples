   
premask <-
function(slicedata, subsamp=TRUE)
{
    niislice <- slicedata$niislice
    mask <- slicedata$mask
    if(subsamp) { # subsample by a factor of 2
        nr <- nrow(niislice)   
        nc <- ncol(niislice)   
        niislice <- niislice[seq(1,nr,2), seq(1,nc,2)]
        mask <- mask[seq(1,nr,2), seq(1,nc,2)]
        nrow <- nrow(niislice)   
        ncol <- ncol(niislice)   
     }
     ## cat("image dimension for simulation:", nrow, ncol, "\n")
     nonactive.level <- mean(range(niislice))
     niislice <- as.vector(niislice) + nonactive.level # required before applying mask
     mask <- as.vector(mask)
     xx <- mask*niislice 
     xv <- xx[xx != 0] # take-off background
     invisible(xv)
}

