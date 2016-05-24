#' Inside region in lumpsucker fishery
#' 
#' Finds regulation area from a position given in dataframe
#' 
#' to be added later
#' 
#' @param data Data frame with columns \code{lat} and \code{lon} givin postions
#' @return Original data with column area added, with index to list of
#' regulation areas in \code{\link{reg.lump}}.
#' @note Could be repeated for other area divisions.
#' @author STJ
#' @seealso \code{\link[geo]{inside.reg.bc}}, \code{\link{reg.lump}}.
#' @references Borrows from \code{\link[geo]{inside.reg.bc}}.
#' @keywords manip
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function (data) 
#' {
#'     if (nrow(data) > 1) 
#'         tmpdata <- data[, c("lat", "lon")]
#'     else tmpdata <- as.data.frame(data[, c("lat", "lon")])
#'     tmpdata$area <- rep(0, nrow(tmpdata))
#'     i <- 1
#'     ind <- geoinside(tmpdata, reg = reg.lump[[i]], option = 0, 
#'         robust = FALSE)
#'     if (length(ind) > 0) 
#'         tmpdata[ind, "area"] <- i
#'     i <- 2
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 3
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 4
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 5
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 6
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 7
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     i <- 8
#'     j <- tmpdata$area == 0
#'     j1 <- c(1:length(j))
#'     j1 <- j1[j == T]
#'     if (length(j1) > 0) {
#'         ind <- geoinside(tmpdata[j1, ], reg = reg.lump[[i]], 
#'             option = 0, robust = FALSE)
#'         if (length(ind) > 0) 
#'             tmpdata[j1[ind], "area"] <- i
#'     }
#'     data$area <- tmpdata$area
#'     return(data)
#'   }
#' 
#' @export inside.reg.lump
inside.reg.lump <-
function (data) 
{
    if (nrow(data) > 1)
        tmpdata <- data[, c("lat", "lon")]
    else tmpdata <- as.data.frame(data[, c("lat", "lon")])
    tmpdata$area <- rep(0, nrow(tmpdata))
    i <- 1
    ind <- geoinside(tmpdata, reg = geo::reg.lump[[i]], option = 0, 
        robust = FALSE)
    if (length(ind) > 0) 
        tmpdata[ind, "area"] <- i
    i <- 2
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 3
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 4
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 5
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 6
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 7
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    i <- 8
    j <- tmpdata$area == 0
    j1 <- c(1:length(j))
    j1 <- j1[j == T]
    if (length(j1) > 0) {
        ind <- geoinside(tmpdata[j1, ], reg = geo::reg.lump[[i]], 
            option = 0, robust = FALSE)
        if (length(ind) > 0) 
            tmpdata[j1[ind], "area"] <- i
    }
    data$area <- tmpdata$area
    return(data)
}
