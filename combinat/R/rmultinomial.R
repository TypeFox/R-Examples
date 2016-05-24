"rmultinomial"<-
function(n, p, rows = max(c(length(n), nrow(p))))
{
# 19 Feb 1997 (John Wallace, 17 Feb 1997 S-news)
# Generate random samples from multinomial distributions, where both n
# and p may vary among distributions
#
# Modified by Scott Chasalow
#
	rmultinomial.1 <- function(n, p)
	{
		k <- length(p)
		tabulate(sample(k, n, replace = TRUE, prob = p), nbins = k)
	}
	#assign("rmultinomial.1", rmultinomial.1)#, frame = 1)
	n <- rep(n, length = rows)
	p <- p[rep(1:nrow(p), length = rows),  , drop = FALSE]
	#assign("n", n)#, frame = 1)
	#assign("p", p)#, frame = 1)
	t(apply(matrix(1:rows, ncol = 1), 1, function(i)
	rmultinomial.1(n[i], p[i,  ])))
}

