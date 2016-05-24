age.window <- function(dat, window,
                       surv = c("enter", "exit", "event")){
    
    if (!is.data.frame(dat))stop("dat must be a data frame")
    if (length(surv) != 3) stop("surv must have length 3")
    fixed.names <- names(dat)
    surv.indices <- match(surv, fixed.names)
    if (length(which(is.na(surv.indices)))){
        x <- which(is.na(surv.indices))
        stop(paste(surv[x], " is not a name in the data frame."))
    }
    
    enter <- dat[, surv.indices[1]]
    exit <- dat[, surv.indices[2]]
    event <- dat[, surv.indices[3]]
    
    who <- (exit > window[1]) & (enter < window[2])
    if (sum(who) > 0.5){ # Not empty selection!
        enter <- enter[who]
        exit <- exit[who]
        event <- event[who]
        ##event <- ifelse(exit > window[2], 0, event)
        event[exit > window[2]] <- 0
        ##exit <- pmin(exit, window[2])
        exit[exit > window[2]] <- window[2]
        ##enter <- pmax(enter[who], window[1])
        enter[enter < window[1]] <- window[1]
        
        dat <- dat[who, ]
        dat[, surv.indices[1]] <- enter
        dat[, surv.indices[2]] <- exit
        dat[, surv.indices[3]] <- event
    }else{ # Empty selection...
        warning(paste("The age interval", window[1], "-", window[2],
                      "is empty."))
        dat <- NULL
    }
    dat
}
