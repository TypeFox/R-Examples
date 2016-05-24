normal.qu <-
function (n = 15, lower = -4, upper = 4, mu = 0, sigma = 1) 
	{
	    qp = seq(lower, upper, length.out = n)
	        qw = dnorm(qp, 0, 1)
	        qw = qw/sum(qw)
	        qp = qp * sigma + mu
	    	
	    return(list(quad.points = qp, quad.weights = qw))
	}

