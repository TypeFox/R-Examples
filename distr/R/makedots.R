.makedotsLowLevel <- function(dots){
       dots$sub <- dots$xlab <- dots$ylab <- dots$main <- dots$type <- NULL
       dots$xlim <- dots$ylim <- dots$yaxt <- dots$axes <- dots$xaxt <- NULL
       dots$panel.last <- dots$panel.first <- dots$frame.plot <- dots$ann <-NULL
       dots$log <- dots$bg <- dots$asp <- NULL
       return(dots)
}
.makedotsL <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots$pch <- dots$cex <- NULL
    dots
}
.makedotsP <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots$lwd <- dots$lty <- NULL
    dots
}
.makedotsPt <- function(dots){
    dots <- dots[names(dots) %in% c("bg", "lwd", "lty")]
    dots["pch"] <- NULL
    if (length(dots) == 0 ) dots <- NULL
    dots
}
