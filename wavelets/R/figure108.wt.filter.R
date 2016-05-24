figure108.wt.filter <- function(filter.objects, level = 1, l = NULL, wavelet = TRUE) {
    draw.leftlabels <- function(y.llabs, left.usrlabelpos)
    {
        for (i in 1:length(y.llabs)) 
        {
            text(left.usrlabelpos, (length(list.filter)-i)*rangeperplot, y.llabs[[i]], xpd = TRUE, cex = .7)
        }
    }

    if(class(filter.objects) != "list") {
        if(class(filter.objects) == "wt.filter") {
            if(filter.objects@wt.name != "none") {
                filter.objects <- wt.filter(filter.objects@wt.name, level = level)
            }
            else {
                filter.objects <- wt.filter(filter.objects@h, level = level)
            }
        }
        if(class(filter.objects) == "character") {
            filter.objects <- wt.filter(tolower(filter.objects), level = level)  
        } 
        if(class(filter.objects) == "numeric") {
            filter.objects <- wt.filter(filter.objects, level = level)
        }
         
        filter.objects <- list(filter.objects)
    }

    for(i in 1:length(filter.objects)) {

        if(class(filter.objects[[i]]) == "numeric") {
            filter.objects[[i]] <- wt.filter(filter.objects[[i]])
        }
        if(class(filter.objects[[i]]) == "character") {
            filter.objects[[i]] <- wt.filter(tolower(filter.objects[[i]]))               
        }
        if(class(filter.objects[[i]]) != "wt.filter") {
            stop("Invalid argument: Every element of filter.objects must be of class 'numeric','character', or 'wt.filter'.")
        }
    }

    if(wavelet == TRUE) {
        list.filter <- list(filter.objects[[1]]@h)

        if(length(filter.objects) != 1) {
            for(i in 2:length(filter.objects)) {
                list.filter <- c(list.filter, list(filter.objects[[i]]@h))
            }
        }      
    }   
    else {
        list.filter <- list(filter.objects[[1]]@g)
   
        if(length(filter.objects) != 1) {
            for(i in 2:length(filter.objects)) {
                list.filter <- c(list.filter, list(filter.objects[[i]]@g))
            }
        }      
    } 

    if(is.null(l) == TRUE) {
        l <- 0:(length(list.filter[[1]]) - 1)
        if(length(list.filter) != 1) {
            for(i in 2:length(list.filter)) {
                if(max(l) < length(list.filter[[i]])) { 
                    l <- 0:(length(list.filter[[i]]) - 1)
                }
            }
        }      
    }
    else {
        if(class(l) == "numeric" && length(l) > 1) {
            l <- max(l)
        }

        l <- 0:l
        for(i in 1:length(list.filter)) {
            if(max(l) < length(list.filter[[i]])) { 
                l <- 0:(length(list.filter[[i]]) - 1)
            }
        }
    } 

    for(i in 1:length(list.filter)) {       
        if(max(l) > (length(list.filter[[i]]) - 1)) {
            add.NA <- rep(NA, times = max(l) - (length(list.filter[[i]]) - 1))
            list.filter[[i]] <- c(list.filter[[i]], add.NA)
        }
    }

    min.list <- 0
    max.list <- 0
  
    for(i in 1:length(list.filter)) {
        if(min.list > min(list.filter[[i]], na.rm = TRUE)) {
             min.list <- min(list.filter[[i]], na.rm = TRUE)
        }        
    }

    for(i in 1:length(list.filter)) {
        if(max.list < max(list.filter[[i]], na.rm = TRUE)) {
             max.list <- max(list.filter[[i]], na.rm = TRUE)
        }        
    }

    plotarray <- 1:length(levels)
    if(abs(max.list) >= abs(min.list)) {
        rangeperplot <- 2*max.list
    }
    else {
        rangeperplot <- 2*min.list
    }
    jstartline <- min.list
    totalrange <- length(list.filter)*(rangeperplot)

    for(j in plotarray) {
        jstartline <- jstartline + rangeperplot
        plotarray[j] <- jstartline
    }

    plot(1, 1, type = "n", xlim = c(0, max(l)), ylim = c(min.list, totalrange + min.list), axes = FALSE, frame.plot = TRUE, xlab = "", ylab = "")
    axis(1) 
   
    for(i in 1:length(list.filter)) {
        par(new = TRUE) 
        shifty1 <- (length(list.filter)-i)*rangeperplot                    
        shifty2 <- (length(list.filter)-i)*rangeperplot + list.filter[[i]] 
        abline(h = shifty1, lty = 1)
        segments(l, shifty1, l, shifty2)
        points(l, shifty2, pch = 15)
    }
 
   y.llabs <- rep(NA, times = length(filter.objects))
   for(i in 1:length(filter.objects)) {
       if(filter.objects[[i]]@wt.name == "haar") {
           y.llabs[i] <- "Haar" 
       }
       else {
           y.llabs[i] <- toupper(filter.objects[[i]]@wt.name)
       }
   }

   draw.leftlabels(y.llabs, par()$usr[1] - abs(par()$usr[1]))

}
