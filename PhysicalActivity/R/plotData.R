plotData <-
function(data, day=NULL, start=NULL, end=NULL){
	findMidnight <- function(data){
	n <- length(data[,1])
		mm <- 0
		for(i in 1:(n-1)){
			mm <- c(mm, ifelse(data$days[i]==data$days[i+1], 0, 1))
			}
		data.midnight <- data[mm==1,]
		first.rowname <- as.numeric(row.names(data[1,]))
		midnightStart <- as.numeric(row.names(data.midnight) ) - first.rowname + 1
		return(midnightStart)
	}
	if(is.null(start)==1){
		start <- 1
		}
	if(is.null(end)==1){
		end <- length(data[,1])
		}
	if(is.null(day)!=1){
		dd <- data[data$days==day,]
		midnightMark <- findMidnight(dd)
		} else{
			dd <- data[start:end,]
			midnightMark <- findMidnight(dd)
			}
	plot(dd$counts, type="l", xlab="Time", ylab="Counts")
	abline(v=midnightMark, lty=2, lwd=1.5, col=4)
	text(midnightMark, 0, pos=1,"0 AM", cex=0.8, col=4)
	}

