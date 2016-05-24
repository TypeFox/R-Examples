bucket.calc <- function(sig, times) {
	start <- which(sig == 1)[1]
	if(is.na(start)) return(NULL)
	start.time <- times[start]
	
	sig <- sig[(start + 1):length(sig)]
	times <- times[(start + 1):length(times)]
	
	if(length(sig) == 0) return(matrix(c(start.time, start.time), ncol = 2))
	
	end <- which(sig == 0)[1]
	if(is.na(end)) return(matrix(c(start.time, tail(times, 1)), ncol = 2))
	end.time <- times[end - 1]
	
	sig <- sig[(end):length(sig)]
	times <- times[(end):length(times)]
	
	rbind(c(start.time, end.time), bucket.calc(sig, times))
}