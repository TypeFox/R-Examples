month.stats.int <-
function(data, fun, ts, years, digits) {
### internal funtion for calculation of month.stats

	m.stats <- matrix(nrow=14, ncol=length(years)+2)
	
	if(fun=="mean") {
		for(y in 1:length(years)) {
			for(m in 1:12) {
				m.stats[m,y] <- mean(data[ts$year==years[y]-1900 & ts$mon==m-1], na.rm=TRUE)
				m.stats[m,length(years)+1] <- mean(data[ts$mon==m-1], na.rm=TRUE)
				m.stats[m,length(years)+2] <- mean(m.stats[m,1:length(years)], na.rm=TRUE)
			}
			m.stats[13,y] <- mean(data[ts$year==years[y]-1900], na.rm=TRUE)
			m.stats[14,y] <- mean(m.stats[1:12,y], na.rm=TRUE)
		}
		m.stats[13,length(years)+1] <- mean(data, na.rm=TRUE)
		m.stats[14,length(years)+2] <- mean(m.stats[1:12,length(years)+2], na.rm=TRUE)
	}
	
	if(fun=="median") {
		for(y in 1:length(years)) {
			for(m in 1:12) {
				m.stats[m,y] <- median(data[ts$year==years[y]-1900 & ts$mon==m-1], na.rm=TRUE)
				m.stats[m,length(years)+1] <- median(data[ts$mon==m-1], na.rm=TRUE)
				m.stats[m,length(years)+2] <- median(m.stats[m,1:length(years)], na.rm=TRUE)
			}
			m.stats[13,y] <- median(data[ts$year==years[y]-1900], na.rm=TRUE)
			m.stats[14,y] <- median(m.stats[1:12,y], na.rm=TRUE)
		}
		m.stats[13,length(years)+1] <- median(data, na.rm=TRUE)
		m.stats[14,length(years)+2] <- median(m.stats[1:12,length(years)+2], na.rm=TRUE)
	}
	
	if(fun=="min") {
		for(y in 1:length(years)) {
			for(m in 1:12) {
				m.stats[m,y] <- m.stats[m,length(years)+1] <- m.stats[m,length(years)+2] <- NA
				if(length(data[ts$year==years[y]-1900 & ts$mon==m-1])>0) m.stats[m,y] <- min(data[ts$year==years[y]-1900 & ts$mon==m-1], na.rm=TRUE)
				if(length(data[ts$mon==m-1])>0) m.stats[m,length(years)+1] <- min(data[ts$mon==m-1], na.rm=TRUE)
				if(any(!is.na(m.stats[m,1:length(years)]))) m.stats[m,length(years)+2] <- min(m.stats[m,1:length(years)], na.rm=TRUE)
			}
			m.stats[13,y] <- NA
			if(length(data[ts$year==years[y]-1900]>0)) m.stats[13,y] <- min(data[ts$year==years[y]-1900], na.rm=TRUE)
		}
		m.stats[13,length(years)+1] <- NA
		if(length(data)>0) m.stats[13,length(years)+1] <- min(data, na.rm=TRUE)
	}
	
	if(fun=="max") {
		for(y in 1:length(years)) {
			for(m in 1:12) {
				m.stats[m,y] <- m.stats[m,length(years)+1] <- m.stats[m,length(years)+2] <- NA
				if(length(data[ts$year==years[y]-1900 & ts$mon==m-1])>0) m.stats[m,y] <- max(data[ts$year==years[y]-1900 & ts$mon==m-1], na.rm=TRUE)
				if(length(data[ts$mon==m-1])>0) m.stats[m,length(years)+1] <- max(data[ts$mon==m-1], na.rm=TRUE)
				if(any(!is.na(m.stats[m,1:length(years)]))) m.stats[m,length(years)+2] <- max(m.stats[m,1:length(years)], na.rm=TRUE)
			}
			m.stats[13,y] <- NA
			if(length(data[ts$year==years[y]-1900]>0)) m.stats[13,y] <- max(data[ts$year==years[y]-1900], na.rm=TRUE)
		}
		m.stats[13,length(years)+1] <- NA
		if(length(data)>0) m.stats[13,length(years)+1] <- max(data, na.rm=TRUE)
	}
	
	if(fun=="sd") {
		for(y in 1:length(years)) {
			for(m in 1:12) {
				m.stats[m,y] <- sd(data[ts$year==years[y]-1900 & ts$mon==m-1], na.rm=TRUE)
				m.stats[m,length(years)+1] <- sd(data[ts$mon==m-1], na.rm=TRUE)
			}
			m.stats[13,y] <- sd(data[ts$year==years[y]-1900], na.rm=TRUE)
		}
		m.stats[13,length(years)+1] <- sd(data, na.rm=TRUE)
	}
	
	m.stats.df <- data.frame(m.stats, row.names=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec",fun,paste0(fun, ".of.months")))
	names(m.stats.df) <- c(years, fun, paste0(fun, ".of.months"))
	
	if(fun=="min" || fun=="max" || fun=="sd") m.stats.df <- m.stats.df[1:13,1:(y+1)]

	for(i in 1:length(m.stats.df)) m.stats.df[,i][is.nan(m.stats.df[,i]) | m.stats.df[,i]==0] <- NA
	
	m.stats.df <- round(m.stats.df, digits)
	return(m.stats.df)
}
