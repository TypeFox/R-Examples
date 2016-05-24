###############################################################################
## Quality control of RFLP data
###############################################################################

RFLPqc <- function(x, rm.band1 = TRUE, QC.lo = 0.8, QC.up = 1.07, QC.rm = FALSE){
    if(length(QC.lo) > 1){
        warning("Only first element of 'QC.lo' is used.")
        QC.lo <- QC.lo[1]
    }
    if(length(QC.up) > 1){
        warning("Only first element of 'QC.up' is used.")
        QC.up <- QC.up[1]
    }
    if(QC.lo <= 0 || QC.lo >= 1)
        stop("'QC.lo' has to be in (0,1).")
    if(QC.up <= 1)
        stop("'QC.up' has to be larger than 1.")        
    temp <- split(x$MW, x$Sample)
    if(QC.rm) rm.samp <- NULL
    for(i in seq_along(temp)){
        tempi <- temp[[i]]
        if((sum(tempi[-1]) < QC.lo*tempi[1]) || (sum(tempi[-1]) > QC.up*tempi[1])){
            cat("Sum of bands of sample ", names(temp)[i], " out of range (", 
                round(sum(tempi[-1])/tempi[1]*100, 2), "%)!", sep = "", fill = TRUE)
            if(QC.rm){
                rm.samp <- c(rm.samp, names(temp)[i])
            }
        }
    }
    
    if(rm.band1){ 
        x <- x[x$Band != 1, ]
        x$Band <- x$Band - 1
    }
    
    if(QC.rm) x <- x[!(x$Sample %in% rm.samp), ]
    
    if(rm.band1 || QC.rm) rownames(x) <- 1:nrow(x)
    
    x
}
