"Mstep.lnorm" <-
function (x, cond, pm, pn) 
{
    nms <- sort(names(pm))
    n <- length(x)
    m <- ncol(cond$u)
    if (all(nms == c("meanlog", "sdlog"))) {
        meanlog <- as.numeric(matrix(log(x), nrow = 1) %*% cond$u)/apply(cond$u, 
            MARGIN = 2, FUN = sum)
        sdlog = sqrt(apply((matrix(log(x), nrow = n, ncol = m) - matrix(meanlog, 
            nrow = n, ncol = m, byrow = TRUE))^2 * cond$u, MARGIN = 2, 
            FUN = sum)/apply(cond$u, MARGIN = 2, FUN = sum))
        return(list(meanlog = meanlog, sdlog = sdlog))
    }
    if (all(nms == "meanlog")) {
        meanlog <- as.numeric(matrix(log(x), nrow = 1) %*% cond$u)/apply(cond$u, 
            MARGIN = 2, FUN = sum)
        return(list(meanlog = meanlog))
    }
    if (all(nms == "sdlog")) {
        sdlog = sqrt(apply((matrix(log(x), nrow = n, ncol = m) - matrix(pn$meanlog, 
            nrow = n, ncol = m, byrow = FALSE))^2 * cond$u, MARGIN = 2, 
            FUN = sum)/apply(cond$u, MARGIN = 2, FUN = sum))
        return(list(sdlog = sdlog))
    }
    stop("Invalid specification of parameters")
}
