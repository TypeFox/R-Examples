pcws2_reg <-
function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c in 2:(C-1))
		{
			x1 <- x[1:c]
			y1 <- y[1:c]
			x2 <- x[c:C]
			y2 <- y[c:C]
			
			a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
			b1 <- -a1 * mean(x1) + mean(y1)

			a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
			b2 <- -a2 * mean(x2) + mean(y2)
			
			ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2)

			if (ss < ssBest) 
				{
					ssBest <- ss
					cBest <- c
					a1Best <- a1
					a2Best <- a2
					b1Best <- b1
					b2Best <- b2
				}
		}
	
	return(list(c=cBest, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2)))
}
