reSurv <- function(time1, id, event, time2 = NULL) {
    if (sum(time1 <= 0) > 0 & is.null(time2)) 
        stop("Observation time must be positive.")
    if (!is.numeric(time1)) 
        stop("Time variable is not numeric")
    if (sum(time1 < 0) > 0 & !is.null(time2)) 
        stop("Observation time must be positive.")
    if (sum(is.wholenumber(id)) < length(id))
        stop("ID must be integer values.")
    if (!is.numeric(time2) & !is.null(time2)) 
        stop("Time variable is not numeric")
    if (length(event) != length(time1)) 
        stop("Time and status are different lengths")
    if (is.logical(event)) 
        status <- as.numeric(event)
    else if (is.numeric(event)) {
        temp <- (event == 0 | event == 1)
        event <- ifelse(temp, event, NA)
        if (sum(is.na(event)) > 0) 
            stop("Invalid status value (Status must be 0 or 1)")
    }
    else
        stop("Invalid status value, must be logical or numeric")
    if (is.null(time2)) 
        ## rc <- data.frame(id = id, Time = time1, event = event)
        rc <- cbind(id = id, Time = time1, event = event)
    if (!is.null(time2)) {
        ## rc <- data.frame(id = id, Time = unlist(lapply(split(time2 - time1, id), cumsum)), event = event)
        rc <- cbind(id = id, Time = unlist(lapply(split(time2 - time1, id), cumsum)), event = event)
    }
    class(rc) <- "reSurv"
    invisible(rc)
}

is.reSurv <- function(x) inherits(x, "reSurv")
is.reReg <- function(x) inherits(x, "reReg")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
