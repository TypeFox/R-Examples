Kullback <- function(Y,X)
{
X <- factor(X)
n <- length(X)
# n.j <- as.matrix(summary(X))
n.j <- summary(X)
p <- ncol(Y)
g <- nlevels(X)
lev.X <- levels(X)
#
# Vj will contain the within-group dispersion matrices of individual groups
Vj <- vector("list", g)
log.ratios <- NA
det.V.j <- NA
# Compute |V|
V <- matrix(0,p,p)
for(k in 1:g)  V <- V + cov(Y[X==lev.X[k],])*(n.j[k]-1)
V <- V/(n-g)
det.V <- det(V)
# Compute Kulback chi-square
Chi <- 0
for(k in 1:g) {
	V.group <- cov(Y[X==lev.X[k],])
	det.V.gr <- det(V.group)
	log.rat <- log(det.V/det.V.gr)
	Chi <- Chi + log.rat*((n.j[k]-1)/2)
	Vj[k] <- list(V.group)
	log.ratios <- c(log.ratios, log.rat)
	det.V.gr <- c(det.V.j, det.V.gr)
	}
#
p.val <- pchisq(Chi,(g-1)*p*(p+1)/2,lower.tail=FALSE)
res <- data.frame(Chi=Chi,df=(g-1)*p*(p+1)/2,p.val=p.val)
names(res) <- c("Chi*","df","P(Chi>Chi*)")
head<-"Kullback test for equal covariance matrices"
out <- list(Kullback=res, V=V, Vj=Vj, log.ratios=log.ratios[-1], det.V=det.V, det.V.j=det.V.j[-1], head=head)
class(out)<-"kback"
out
}


print.kback<-function(x,digits= max(3, getOption("digits")), ...){
cat("\n")
cat(x$head,"\n")
rq<-x$Kullback
print(rq,digits=digits)
cat("\n")
invisible(x)
}