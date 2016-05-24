figure98.wt.filter <- function(filter, levels = NULL, wavelet = TRUE, y.normalize = TRUE) {
    draw.leftlabels <- function(y.llabs, left.usrlabelpos)
    {
        for (i in 1:length(y.llabs)) 
        {
            text(left.usrlabelpos, (length(list.filter)-i)*rangeperplot + abs(min.list), y.llabs[[i]], xpd = TRUE)
        }
    }

    if(is.null(levels) == TRUE) {
        levels <- 1:7
    }
    if(class(levels) == "numeric" && length(levels) == 1) {
        levels <- 1:levels
    }

    if(class(filter) == "numeric") {
        if(length(filter)%%2 != 0) {
            stop("Invalid argument: filter length must be even.")
        } 
        filter.name <- ""
    }
    if(class(filter) == "character") {
        if(tolower(filter) != "haar") {
            filter.name <- toupper(filter)
        }
        else {
            filter.name <- "Haar"
        }
        filter <- tolower(filter)
    }
    if(class(filter) == "wt.filter") {
        filter.name <- toupper(filter@wt.name)
        filter <- filter@h
    }

    if(wavelet == TRUE) {
        list.filter <- list(wt.filter(filter, level = levels[1])@h)

        if(length(levels) != 1) {
            for(i in 2:length(levels)) {
                list.filter <- c(list.filter, list(wt.filter(filter, level = levels[i])@h))
            }
        }      
    }   
    else {
        list.filter <- list(wt.filter(filter, level = levels[1])@g)
   
        if(length(levels) != 1) {
            for(i in 2:length(levels)) {
                list.filter <- c(list.filter, list(wt.filter(filter, level = levels[i])@g))
            }
        }      
    } 

    min.list <- 0
    max.list <- 0
   
    for(i in 1:length(list.filter)) {
        if(min.list > min(list.filter[[i]])) {
             min.list <- min(list.filter[[i]])
        }        
    }

    for(i in 1:length(list.filter)) {
        if(max.list < max(list.filter[[i]])) {
             max.list <- max(list.filter[[i]])
        }        
    }

    plotarray <- 1:length(levels)
    rangeperplot <- 2*max.list
    jstartline <- min.list
    totalrange <- length(list.filter)*(rangeperplot)

    if(y.normalize == TRUE) {
        for(i in 1:length(list.filter)) {
            if(max.list - abs(max(list.filter[[i]])) < max.list - abs(min(list.filter[[i]]))) {
                scale <- max.list/abs(max(list.filter[[i]]))
            }
            else {
                scale <- max.list/abs(min(list.filter[[i]]))
            }
          
            list.filter[[i]] <- list.filter[[i]]*scale
        }
    }
   
    for(j in plotarray) {
        jstartline <- jstartline + rangeperplot
        plotarray[j] <- jstartline
    }

    plot(0, 0, type = "n", xlim = c(0, length(list.filter[[1]])-1), ylim = c(min.list, totalrange),  axes = FALSE, frame.plot = TRUE, xlab = "l", ylab = "j")
 
    for(i in 1:length(list.filter)) {
        par(new = TRUE) 
        plot(0, 0, type = "n", xlim = c(0, length(list.filter[[i]])-1), ylim = c(min.list, totalrange), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")            
        par(xaxp = c(0,length(list.filter[[i]])-1,1))
        axis(1, at = axTicks(1), labels = c(0, expression(L[j]-1)))
        shifty <- (length(list.filter)-i)*rangeperplot + abs(min.list)            
        lines(0:(length(list.filter[[i]])-1), shifty + list.filter[[i]], type = "l")
        abline(h = shifty, lty = 2)
    }

    y.llabs <- as.character(levels)   

    draw.leftlabels(y.llabs, par()$usr[1] - abs(par()$usr[1]))
}












