rz.transform <-
function (x, jitter = FALSE){
    
    y <- x[!is.na(x)]
    mean.y <- mean(y, na.rm = TRUE)
    sd.y <- sd(y, na.rm = TRUE)

	#if there are infinite values in y
	#set these to 10 less than/greater than
	#the real min and max of y
    y[y == Inf] <- max(y[y < Inf]) + 10
    y[y == -Inf] <- min(y[y > -Inf]) - 10

    if(jitter){
        y <- rank(y + runif(length(y))/(sd(y) * 10^8))
        }else{
        	y <- rank(y)
        	}
    
    x[!is.na(x)] <- qnorm((y - 0.5)/length(y))
    
    return(x * sd.y/sd(x, na.rm = TRUE) - mean(x, na.rm = TRUE) + mean.y)
	}
