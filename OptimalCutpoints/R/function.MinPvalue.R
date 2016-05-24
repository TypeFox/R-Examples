function.MinPvalue <-
function(data, marker, status, tag.healthy = 0, direction = c("<", ">"), control = control.cutpoints(), pop.prev, ci.fit = FALSE, conf.level = 0.95, measures.acc){
	direction <- match.arg(direction)
	
	marker.healthy = data[data[,status] == tag.healthy, marker]
	marker.diseased = data[data[,status] != tag.healthy, marker]

	n <- length(measures.acc$cutoffs)-1

	stat <- numeric(n)
	pvalue <- numeric(n)
	RR <- numeric(n)

	if (direction == "<") 
	{
		for (i in 1:n) 
		{
			tabl <- matrix(NA,nrow=2,ncol=2)
			tabl[1,1] <- apply(outer(marker.diseased, measures.acc$cutoffs[i+1],">="),2,sum)
			tabl[1,2] <- apply(outer(marker.healthy, measures.acc$cutoffs[i+1],">="),2,sum)
			tabl[2,1] <- apply(outer(marker.diseased, measures.acc$cutoffs[i+1],"<"),2,sum)
			tabl[2,2] <- apply(outer(marker.healthy, measures.acc$cutoffs[i+1],"<"),2,sum)

			test <- chisq.test(tabl)
			stat[i] <- test$statistic
			pvalue[i] <- test$p.value			
 
			if (any(tabl == 0)) tabl = tabl+0.5
			RR[i] <-( (tabl[1,1])/(tabl[1,1]+tabl[1,2])) / ((tabl[2,1])/(tabl[2,1]+ tabl[2,2]))	
		}
	}


	if (direction == "<") 
	{
		for (i in 1:n) 
		{
			tabl <- matrix(NA,nrow=2,ncol=2)
			tabl[1,1] <- apply(outer(marker.diseased, measures.acc$cutoffs[i+1],"<"),2,sum)
			tabl[1,2] <- apply(outer(marker.healthy, measures.acc$cutoffs[i+1],"<"),2,sum)
			tabl[2,1] <- apply(outer(marker.diseased, measures.acc$cutoffs[i+1],">="),2,sum)
			tabl[2,2] <- apply(outer(marker.healthy, measures.acc$cutoffs[i+1],">="),2,sum)

			test <- chisq.test(tabl)
			stat[i] <- test$statistic
			pvalue[i] <- test$p.value
			
 			if (any(tabl == 0)) tabl = tabl+0.5
			RR[i] <-( (tabl[1,1])/(tabl[1,1]+tabl[1,2])) / ((tabl[2,1])/(tabl[2,1]+ tabl[2,2]))	
		}
	}


	# Different methods for adjusting the minimum p-value:												   
	lower <- measures.acc$cutoffs[2]
	upper <- measures.acc$cutoffs[n+1]
	epsi.high <- (length(which(data[,marker]>=upper))/length(data[,marker]))*100
	epsi.low <- (length(which(data[,marker]<=lower))/length(data[,marker]))*100

	# Miller and Siegmund's formula for adjusting the minimum p-value:
	PADJMS <- function(cutpoint,pvalue,epsi.high,epsi.low) {
		pmin <- min(pvalue, na.rm = TRUE)		
		cut.point <- cutpoint[which(round(pvalue,10) == round(pmin,10))]
	 	z <- qnorm(1-pmin/2)
		f.z <- dnorm(z)
		pacor <- f.z*(z-1/z)*log((epsi.high*(1-epsi.low))/((1-epsi.high)*epsi.low))+(4*f.z)/z
		pval <- c(cut.point,pmin,epsi.high,epsi.low,pacor)
		names(pval) <- c("cutpoint","pmin","epsi.high","epsi.low","pms")
		return(pval)
	}

	# Altman's formula for adjusting the minimum p-value
	# (epsi.high=epsi.low=5%, epsi.high=epsi.low=10%):  
	PALT510 <- function(cutpoint,pvalue) {
		pmin <- min(pvalue, na.rm = TRUE)
		cut.point <- cutpoint[round(pvalue,10) == round(pmin,10)]
		pcor10 <-(-1.63*pmin*(1+2.35*log(pmin)))
		pcor5 <-(-3.13*pmin*(1+1.65*log(pmin)))
		pval <- c(cut.point,pmin,pcor5,pcor10)
		names(pval) <- c("cutpoint","pmin","palt5","palt10")
		return(pval)
	}

	cutpoints <- measures.acc$cutoffs[2:length(measures.acc$cutoffs)]

	if (control$adjusted.pvalue == "PADJMS") {
		cMinPvalue <- PADJMS(cutpoints,pvalue,epsi.high,epsi.low)[1]
		minimum.pvalue <- PADJMS(cutpoints,pvalue,epsi.high,epsi.low)[2]
		minimum.adjusted.pvalue <- PADJMS(cutpoints,pvalue,epsi.high,epsi.low)[5]
	}

	if (control$adjusted.pvalue == "PALT5") {
		cMinPvalue <- PALT510(cutpoints,pvalue)[1]
		minimum.pvalue <- PALT510(cutpoints,pvalue)[2]
		minimum.adjusted.pvalue5 <- PALT510(cutpoints,pvalue)[3]
	}

	if (control$adjusted.pvalue == "PALT10") {
		cMinPvalue <- PALT510(cutpoints,pvalue)[1]
		minimum.pvalue <- PALT510(cutpoints,pvalue)[2]
		minimum.adjusted.pvalue10 <- PALT510(cutpoints,pvalue)[4]
	}	

	optimal.cutoff <- obtain.optimal.measures(cMinPvalue, measures.acc)

	if (control$adjusted.pvalue == "PADJMS") {
		if (is.na(minimum.adjusted.pvalue)) optimal.criterion = minimum.pvalue
		else optimal.criterion = minimum.adjusted.pvalue
		res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = c(NA,pvalue), stat = c(NA,stat), RR = c(NA,RR), minimum.pvalue = minimum.pvalue, optimal.criterion = optimal.criterion)
	}

	if (control$adjusted.pvalue == "PALT5") {
		if (is.na(minimum.adjusted.pvalue5)) optimal.criterion = minimum.pvalue
		else optimal.criterion = minimum.adjusted.pvalue5
		res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = c(NA,pvalue), stat = c(NA,stat), RR = c(NA,RR), minimum.pvalue = minimum.pvalue, optimal.criterion = optimal.criterion)
	}

	if (control$adjusted.pvalue == "PALT10") {
		if (is.na(minimum.adjusted.pvalue10)) optimal.criterion = minimum.pvalue
		else optimal.criterion = minimum.adjusted.pvalue10
		res <- list(measures.acc = measures.acc, optimal.cutoff = optimal.cutoff, criterion = c(NA,pvalue), stat = c(NA,stat), RR = c(NA,RR), minimum.pvalue = minimum.pvalue, optimal.criterion = optimal.criterion)
	}

	res
}
