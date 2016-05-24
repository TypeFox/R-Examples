plot.emg <- function(x, channels = "all", samples = 0, type = "l", timeunits = c("samples", 
    "seconds"), add = FALSE, ...) {
    object <- x
    args <- list(...)
    namesargs <- names(args)
    timeunits <- match.arg(timeunits)
    if ((timeunits == "seconds") & (object$samplingrate <= 0)) 
        warning("Select timeunits='samples' or provide the sampling rate to determine the total duration of the EMG")
    if (!is.vector(object$values) & !is.vector(channels)) 
        stop("'channels' should be an integer, a string or a vector")
    if (!is.numeric(samples) | !is.vector(samples) | !(length(samples) %in% c(1, 
        2))) 
        stop("Not valid 'samples' parameter")
    if (is.vector(object$values)) {
        x <- 1:length(object$values)
    } else {
        x <- 1:dim(object$values)[1]
    }
    
    if (length(samples) == 1) {
        if (samples != 0) {
            if (samples < 0) 
                stop("samples (Number of samples) need to be a positive integer value")
            if (samples > length(object$values)) {
                warning("Required number of samples to be plot is larger than the length of data")
            } else {
                x <- head(x, samples)
                object$values <- head(object$values, samples)
            }
        }
    } else {
        if ((samples[2] < samples[1]) | any(samples < 0) | any(samples > length(object$values))) 
            stop("Samples (range of samples) is not valid")
        x <- x[samples[1]:samples[2]]
        if (is.vector(object$values)) {
            object$values <- object$values[samples[1]:samples[2]]
        } else {
            object$values <- object$values[samples[1]:samples[2], ]
        }
        
    }
    
    if ((timeunits == "seconds") & (object$samplingrate > 0)) {
        x <- x/object$samplingrate
        xlab = "time (seconds)"
    } else xlab = "time (samples)"
    
    if (!is.vector(object$values)) {
        nch <- dim(object$values)[2]
    } else {
        nch <- 1
    }
    
    if (length(channels) == 1 & channels[1] == "all") 
        channels <- 1:nch
    
    if (is.vector(object$values) | length(channels) == 1) {
        if (!is.vector(object$values)) {
            if (is.numeric(channels)) {
                ch <- channels[1]
            } else {
                if (is.character(channels)) {
                  ch <- grep(channels, object$data.name)
                  if (length(ch) == 0) 
                    stop(paste("'", channels, "' is not  a valid channel name"))
                } else {
                  stop(paste(channels, "is not  a valid channel in data"))
                }
            }
            
            if (ch < 1 | ch > dim(object$values)[2]) 
                stop(paste(channels, "is not  a valid channel in data"))
            
            object$values <- object$values[, ch]
            object$data.name <- object$data.name[ch]
            object$units <- object$units[ch]
        }
        if (!("ylab" %in% namesargs)) {
            if (object$data.name != "") 
                ylab <- object$data.name else ylab <- "EMG"
            if (object$units != "") 
                ylab <- paste(ylab, "(", object$units, ")", sep = "")
            if (add) 
                lines(x, object$values, type = type, ylab = ylab, ...) else plot(x, object$values, type = type, ylab = ylab, xlab = xlab, ...)
        } else {
            if (add) 
                lines(x, object$values, type = type, ...) else plot(x, object$values, type = type, xlab = xlab, ...)
        }
    } else {
        if (add) 
            warning("'add' parameter ommited, add to an already existing plot not implemented for multichannel data")
        if (is.numeric(channels)) {
            if (any(channels < 1) | any(channels > nch)) 
                stop(paste("No valid channel numbers, must be between 1 and", nch))
            object$values <- object$values[, channels]
            object$data.name <- object$data.name[channels]
            object$units <- object$units[channels]
        } else {
            if (is.character(channels)) {
                chs <- numeric()
                for (i in 1:nch) {
                  ch <- grep(channels[i], object$data.name)
                  if (length(ch) == 0) 
                    stop(paste("'", channels[i], "' is not  a valid channel name"))
                  chs[i] <- ch
                }
                object$values <- object$values[, chs]
                object$data.name <- object$data.name[chs]
                object$units <- object$units[chs]
            } else {
                stop("No valid channels parameter")
            }
        }
        
        nch <- dim(object$values)[2]
        if ("ylab" %in% namesargs) 
            warning("Parameter 'ylab' not used in multichannel data")
        op <- par(no.readonly = TRUE)
        layout(matrix(1:nch, ncol = 1), widths = 1, heights = c(rep(5, nch - 1), 
            7), respect = FALSE)
        par(mar = c(0, 4, 2, 1), bty = "l")
        for (i in 1:(nch - 1)) {
            ylab <- object$data.name[i]
            if (object$units[i] != "") 
                ylab <- paste(ylab, "(", object$units[i], ")", sep = "")
            if ("main" %in% namesargs & i > 1) {
                args2 <- args
                args2$main <- NULL
                do.call("plot", c(list(x = x, y = object$values[, i], type = type, 
                  ylab = ylab, xaxt = "n"), args2))
                
            } else {
                plot(x, object$values[, i], type = type, ylab = ylab, xaxt = "n", 
                  ...)
            }
        }
        par(mar = c(4, 4, 2, 1))
        ylab <- object$data.name[nch]
        if (object$units[nch] != "") 
            ylab <- paste(ylab, "(", object$units[nch], ")", sep = "")
        if ("main" %in% namesargs) {
            args2 <- args
            args2$main <- NULL
            do.call("plot", c(list(x = x, y = object$values[, nch], type = type, 
                ylab = ylab, xlab = xlab), args2))
            
        } else {
            plot(x, object$values[, nch], type = type, ylab = ylab, xlab = xlab, 
                ...)
        }
        
        op$fig <- NULL
        par(op)
    }
} 
