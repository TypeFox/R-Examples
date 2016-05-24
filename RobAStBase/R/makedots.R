## dots modifications
.makedotsLowLevel <- function(dots){
       dots$sub <- dots$xlab <- dots$ylab <- dots$main <- dots$type <- NULL
       dots$xlim <- dots$ylim <- dots$yaxt <- dots$axes <- dots$xaxt <- NULL
       dots$panel.last <- dots$panel.first <- dots$frame.plot <- dots$ann <-NULL
       dots$log <- dots$asp <- NULL
       return(dots)
}
.deleteDotsABLINE <- function(dots){
    dots$reg <- dots$a <- dots$b <- NULL
    dots$untf <- dots$h <- dots$v <- NULL
    dots
}
.deleteDotsTEXT <- function(dots){
   dots$labels <- dots$offset <- dots$vfont <- dots$pos <- dots$font <- NULL
   dots
}
.makedotsL <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots$pch <- dots$cex <- NULL
    .deleteDotsABLINE(.deleteDotsTEXT(dots))
}
.makedotsP <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots$lwd <- NULL
    .deleteDotsABLINE(.deleteDotsTEXT(dots))
}
.makedotsPt <- function(dots){
      dots <- dots[names(dots) %in% c("bg", "lwd", "lty")]
      if (length(dots) == 0 ) dots <- NULL
      return(dots)
}
.makedotsAB <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots <- .deleteDotsTEXT(dots)
    dots$pch <- dots$cex <- NULL
}
.makedotsT <- function(dots){
    dots <- .makedotsLowLevel(dots)
    dots <- .deleteDotsABLINE(dots)
    dots
}
