mystat <-
function(x, numdig=3, na.rm=TRUE,printit=TRUE)
{
	x<-as.vector(x)
	if (na.rm)   x <- x[!is.na(x)]
	skewx<-theskew(x)
	kurtosisx<-thekurt(x)
	allstats<-data.frame(cbind(min=min(x), max=max(x), mean=mean(x), median=median(x),sdev=sd(x), skew=theskew(x), kurtosis=thekurt(x)),row.names='')
	#too long for one line
	if(printit)
		{
		print(format(allstats[1:4],digits=numdig))
		print(format(allstats[5:7],digits=numdig))
	}
	return(invisible(allstats))
	}
