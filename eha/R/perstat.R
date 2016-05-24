perstat <- function(surv, period, age = c(0, 200)){
    if (ncol(surv) != 4) stop("Need a full 'surv' object with four columns.")
    surv <- data.frame(surv)
    names(surv) <- c("enter", "exit", "event", "birthdate")
    n.rows <- length(period) - 1
    n.cols <- length(age) - 1
    row.name <- character(n.rows)
    for (i in 1:n.rows){
        row.name[i] <-
            paste("(", as.character(period[i]),
                  " - ", as.character(period[i+1]), "]", sep = "")  
    }
    ##if (n.cols > 1){
        col.name <- character(n.cols)
        for (i in 1:n.cols){
            col.name[i] <-
                paste("(", as.character(age[i]),
                      " - ", as.character(age[i+1]), "]", sep = "")  
        }
    ##}

    events <- matrix(0, ncol = n.cols, nrow = n.rows)
    exposure <- matrix(0, ncol = n.cols, nrow = n.rows)
    intensity <- matrix(0, ncol = n.cols, nrow = n.rows)
    rownames(events) <- row.name
    rownames(exposure) <- row.name
    rownames(intensity) <- row.name
    ##if(n.cols > 1){
        colnames(events) <- col.name
        colnames(exposure) <- col.name
        colnames(intensity) <- col.name
    ##}
    
    for (i in 1:n.rows){
        per.dat <- cal.window(surv, c(period[i], period[i + 1]))
        ##if (nrow(per.dat) > 0){
        if (!is.null(per.dat)){
            for (j in 1:n.cols){
                pa.dat <- age.window(per.dat, c(age[j], age[j + 1]))
                ##nr <- nrow(pa.dat)
                if (!is.null(pa.dat)){
                    events[i, j] <- sum(pa.dat$event)
                    exposure[i, j] <- sum(pa.dat$exit - pa.dat$enter)
                    intensity[i, j] <- events[i, j] / exposure[i, j]
                }else{
                    intensity[i, j] <- NaN
                }
            }
        }else{
            for (j in 1:n.cols){
                intensity[i, j] <- NaN
            }
        }
    }
    
    list(events = events,
         exposure = exposure,
         intensity = intensity)
}
          
