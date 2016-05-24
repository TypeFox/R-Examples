pcws3_reg <-
function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c1 in 2:(C-2))
		{
			for (c2 in (c1+1):(C-1))
				{
			
					x1 <- x[1:c1]
					y1 <- y[1:c1]
					x2 <- x[c1:c2]
					y2 <- y[c1:c2] 
					x3 <- x[c2:C]
					y3 <- y[c2:C] 
			
					a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
					b1 <- -a1 * mean(x1) + mean(y1)

					a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
					b2 <- -a2 * mean(x2) + mean(y2)
			
					a3 <- sum((x3-mean(x3))*(y3-mean(y3)))/sum((x3-mean(x3))^2)
					b3 <- -a3 * mean(x3) + mean(y3)

			ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2) + sum((a3*x3+b3-y3)^2)

			if (ss < ssBest) 
				{
					ssBest <- ss
					c1Best <- c1
					c2Best <- c2
					a1Best <- a1
					b1Best <- b1
					a2Best <- a2
					b2Best <- b2
					a3Best <- a3
					b3Best <- b3
				}
		}
	}
	return(list(c1=c1Best, c2=c2Best, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, a3=a3Best, b3=b3Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2,a3*x3+b3-y3)))
}
