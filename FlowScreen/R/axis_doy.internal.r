#' Create custom axis starting on hyrologic year start month
#' 
#' @param hyrstart numeric indicating month for start of the hydrologic year 
#'   (water year).
#' @author Paul Whitfield

axis_doy.internal <- function(hyrstart=10) {
    
    cday <-c(1,32,60,91,121,152,182,213,244,274,305, 335,366,397,425,456,486,517,547,578,609,639,670)
    ctxt <-c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov","Dec",
             "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    
    wday<-cday[hyrstart:(hyrstart+11)]-cday[hyrstart]+1
    wtxt<-ctxt[hyrstart:(hyrstart+11)]
    
    graphics::axis(side=1, at=wday , labels=wtxt, line = 0, tck = -0.025, xlab="",  xlab="")
}