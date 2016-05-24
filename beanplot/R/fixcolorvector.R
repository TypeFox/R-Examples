`fixcolorvector` <-
function(col) {
    mrgb <- function(col) {
        rgb(col[1], col[2], col[3], maxColorValue = 255)
    }
    #mrgba <- function(col, alpha) {
    #    rgb(col[1], col[2], col[3], alpha * 255, maxColorValue = 255)
    #}
    if (missing(col) || ((length(col) == 1) && all(col2rgb(col, 
        alpha = TRUE) == col2rgb(par("fg"), alpha = TRUE)))) {
        col <- c(par("fg"), mrgb(col2rgb(par("bg"))), par("fg"))
        #if (names(dev.cur()) %in% c("pdf", "quartz")) {
        #    col <- c(par("fg"), mrgba(col2rgb(par("bg")), 0.8), 
        #        mrgba(col2rgb(par("fg")), 0.8))
        #}
    }
    else if (length(col) == 1) {
        col <- c(col, par("fg"), par("fg"))
        #if (names(dev.cur()) %in% c("pdf", "quartz")) {
        #    col <- c(col, mrgba(col2rgb(par("fg")), 0.8), mrgba(col2rgb(par("fg")), 
        #        0.8))
        #}
    }
    else if (length(col) == 2) {
        col <- c(col[1], col[2], par("fg"))
        #if (names(dev.cur()) %in% c("pdf", "quartz")) {
        #    col <- c(col[1], col[2], mrgba(col2rgb(par("fg")), 
        #        0.8))
        #}
    }
    if (length(col) == 3) {
        col[4] <- par("fg")
    }
    col
}

