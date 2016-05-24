#' @name segment.print
#' @title Print Segmented Statistics
#' @param mod  Segmented model object
#' @keywords internal

segment.print <- function(mod) {

    LCL <- round(mod$psi[[2]]-(mod$psi[[3]]*1.644854), digits = 4)
    UCL <- round(mod$psi[[2]]+(mod$psi[[3]]*1.644854), digits = 4)

    segout <- matrix(nrow=2, ncol=7)
    segout[1,1] <- "Breakpoint Dose"
    segout[1,3] <- "Lower CL"
    segout[1,5] <- "Upper CL"
    segout[1,7] <- "Std Error"
    
    segout[2,1] <- round(mod$psi[[2]], digits = 4)
    segout[2,3] <- LCL
    segout[2,5] <- UCL
    segout[2,7] <- round(mod$psi[[3]], digits = 4)
    
    drsmooth.print(segout)
    
}