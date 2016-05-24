check.dist <- function(sp, pp, main = NULL, col = NULL, printLegend = TRUE){
    if (!inherits(sp, "coxreg")){
        if (inherits(pp, "coxreg")){ # swap:
            tmp <- pp
            pp <- sp
            sp <- tmp
            rm(tmp)
        }else{
            stop ("Some argument must be of type 'coxreg'")
        }
    }
    if (!inherits(pp, "phreg"))
        stop ("Some argument must be of type 'phreg' or 'pchreg'")

    if (!sp$nullModel){
        if ((!sp$center) && pp$center)
            warning("The non-parametric fit is not centered.") 
        if ((!pp$center) && sp$center)
            warning("The parametric fit is not centered.")
    }
    if ((!is.null(sp$strata)) || (!is.null(pp$strata)))
        stop("Not for stratified fits; try a comparison stratum by stratum.") 
    if (is.null(main)){
        main <- pp$dist # Capitalize:
        substr(main, 1, 1) <- toupper(substr(main, 1, 1))
        if (main == "Pch") main <- "Piecewise constant"
        if (main == "Ev") main = "Extreme value"
    }

    x.max <- max(pp$y[, 2])
    x <- plot.coxreg(sp, fn = "cum", fig = FALSE)
    if (is.null(x)){
        cat("Error: Must be fixed in check.dist!")
        return(x)
    }
    
    x[[1]][, 2] <- cumsum(x[[1]][, 2]) # Added in 2.4-2
    y.max <- max(x[[1]][, 2])
    if (length(x) > 1){
        for (i in 2:length(x)) y.max <- max(y.max, x[[i]][, 2])
    }
    if (is.null(col)){
        col <- c(1, 1)
    }else{
        if (length(col) != 2) stop("Length of 'col' must be 0 or 2.")
    }
    plot(pp, fn = "cum", fig = TRUE, ## Removed 2.4-0: new.data = pp$means,
         ylim = c(0, y.max), main = main, col = col[1])
    for (rr in 1:length(x)){
        xx <- x[[rr]]
        xx <- rbind(xx, c(x.max, xx[NROW(xx), 2])) # Added 2011-08-10 (2.0-3)
        lines(xx[, 1], xx[, 2], type = "s", lty = 2, col = col[2])
    }
    if (printLegend){
        legend(x = "topleft", legend = c("Parametric", "Non-parametric"),
               lty = 1:2, col = col)
    }
}

