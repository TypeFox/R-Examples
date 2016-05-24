summary.rmas <-
function (object, stage = NULL, harvest=FALSE, ...) 
{
   cosa <- object
   if(length(names(cosa))>0){ # bifurcación provisional para separar rmas de projectn y de projectn2
       if(harvest==FALSE) cosa<- cosa$vn else cosa <-cosa$harvest # seleciona trayectorai de la población o el manejo
    }
    nl <- length(cosa)
    time <- dim(cosa[[1]])[2]
    if (nl > 1) {
        abundances <- sapply(cosa, function(x) apply(x, 2, sum))
        if (!is.null(stage)) 
            abundances <- sapply(cosa, function(x) x[stage, ])
        summary <- apply(abundances, 1, function(x) c(min(x), 
            mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x)))
        summary <- t(summary)
        summary <- data.frame(cbind(0:(time - 1), round(summary, 
            2)), row.names = 0:(time - 1))
        names(summary) <- c("Time", "Minimum", "-1 S.D.", "Average", 
            "+1 S.D.", "Maximum")
    }
    if (nl == 1) {
        abundances <- apply(cosa[[1]], 2, sum)
        if (!is.null(stage)) 
            abundances <- apply(cosa, function(x) x[stage, ])
        summary <- cbind(0:(time - 1), abundances)
        colnames(summary) <- c("Time", "Abundance")
    }
    summary <- data.frame(summary, row.names = NULL, check.names = FALSE)
    class(summary) <- c("summary.rmas", class(summary))
    
     warnold <- options("warn")
     options(warn=-1)
         plot(summary)
    options(warn=warnold$warn)
    return(summary)
}

