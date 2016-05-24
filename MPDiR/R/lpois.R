lnorm <- function(p, d) {
	-sum(d$nyes * pnorm((log(d$Q) - p[1])/p[2], 
		lower.tail = TRUE,  log.p = TRUE) + 
	d$nno * pnorm((log(d$Q) - p[1])/p[2], 
		lower.tail = FALSE, log.p = TRUE))
	}
	
lpois <- function(p, d) {
	-sum(d$nyes * ppois(p[2], d$Q/p[1], lower.tail = FALSE, 
		log.p = TRUE) + 
	d$nno * ppois(p[2], d$Q/p[1], lower.tail = TRUE, 
		log.p = TRUE))
	}
	
lpois1 <- function(q, p, d) lpois(c(q, p), d)
