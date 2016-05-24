breaks.AP <-
function(posture)
{
	y <- posture
	n <- length(y)
	
	mmm <- length(y)
	one <- y[-mmm]
	two <- y[-1]
	
	# transitions from sed to not
	trans.up <- 
		(one=="0")&(two!="0")
	num.up.AP <- sum(trans.up,na.rm=T)
	if (num.up.AP==0)
		num.up.AP <- NA
	return(num.up.AP=num.up.AP)	
}

