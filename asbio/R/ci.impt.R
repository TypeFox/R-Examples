ci.impt <- function(y1, n1, y2 = NULL, n2 = NULL, avail.known = FALSE, pi.2 = NULL, conf = .95, x100 = TRUE, alternative = "two.sided", bonf = TRUE, wald = FALSE){
pi.1 <- y1/n1
alpha <- 1 - conf
if(bonf == TRUE) alpha <- alpha/length(pi.1)
z.star <- ifelse(alternative == "two.sided", qnorm(1 - (alpha/2)), qnorm(1 - alpha))

if(avail.known == FALSE){
pi.2 <- y2/n2
sigma.hat.I <- sqrt((1 - pi.1)/(n1 * pi.1) + (1 - pi.2)/(n2 * pi.2))
}

if(avail.known == TRUE){
if(any(pi.2 < 0) | any(pi.2 > 1)) stop("pi.2 must be in the interval (0, 1)")
sigma.hat.I <- sqrt((1 - pi.1)/(n1 * pi.1))
if (wald == TRUE) {
sigma.hat.W <- sqrt(pi.2^2 * (pi.1 * (1 - pi.1))/n1)
if(x100 == TRUE){
sigma.hat.W <- sqrt(100^2 * pi.2^2 * (pi.1 * (1 - pi.1))/n1)
}
}
}

if(x100 == TRUE) theta.hat <- (pi.1*100) * pi.2  else  theta.hat <- pi.1 * pi.2

CI.L <- theta.hat * exp(-1 * z.star * sigma.hat.I)
CI.U <- theta.hat * exp(z.star * sigma.hat.I)
if(avail.known==TRUE & wald == TRUE){
CI.L <- theta.hat - z.star * sigma.hat.W
CI.U <- theta.hat + z.star * sigma.hat.W
}

CI.L <- ifelse(CI.L < 0, 0, CI.L)
CI.U <- ifelse(x100 == TRUE & CI.U > 100, 100, CI.U)
CI.U <- ifelse(x100 == FALSE & CI.U > 1, 1, CI.U)
if(alternative =="less") CI.U <- NA
if(alternative == "greater") CI.L <- NA
CI <- data.frame(theta.hat = theta.hat, CI.L = CI.L, CI.U = CI.U)

head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence interval for product of proportions"))
if(length(pi.1) > 1) head<-paste(paste(as.character(conf*100),"%",sep=""),c("Confidence intervals for product of proportions"))
if(bonf == TRUE & length(pi.1) > 1) head <- paste(head, "\nBonferroni-adjusted intervals")
ends<-c("Estimate",paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
res<-list(SE=ifelse(avail.known ==TRUE & wald==TRUE, sigma.hat.W, sigma.hat.I), ci=CI, ends=ends, head=head)
class(res)<-"ci"
res
}
