"dots2" <-
function( x, y, colx='green', coly='blue',
                  lab1 = deparse(substitute(x)), 
                   lab2 = deparse(substitute(y)), ... ){
		
	sx1 <- sort(x)
	sy1 <- unlist(lapply(table(sx1),seq))
	
	sx2 <- sort(y)
	sy2 <- unlist(lapply(table(sx2),seq))
	
	sy1 <- sy1/ (max(sy1,sy2)+2)
	sy2 <- sy2/ (max(sy1,sy2)+2) + 1
	
	plot( c(sx1,sx2), c(sy1,sy2), xlab="", ylab="", yaxt="n",ylim=c(0,2),
             type="n",...)
	points( sx1, sy1, col=colx,...)
	points( sx2, sy2, col=coly,...)
	axis(2, at=c(0.5,1.5), labels= c(lab1,lab2),srt=90,tick=FALSE)

}

