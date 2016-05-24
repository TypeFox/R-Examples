boot.ci.M<-function(X1, X2, alpha = 0.05, est = huber.one.step, R = 1000){
X1 <- X1[complete.cases(X1)]; X2 <- X2[complete.cases(X2)]

X1b <- bootstrap(X1, R = R, statistic = est)$dist
X2b <- bootstrap(X2, R = R, statistic = est)$dist
bvec <- sort(X1b - X2b)
test <- sum(bvec < 0)/R + 0.5 * sum(bvec == 0)/R
pv = 2 * min(c(test, 1 - test))
low <- round((alpha/2) * R)
up <- round((1 - alpha/2) * R)
se <- sqrt(var(bvec))

head <- paste(paste(as.character((1-alpha)*100),"%",sep=""),c("Confidence interval for true difference of location measures :"))
head2 <- "Test results :"
ends <- c(paste(as.character(c((alpha)/2,1-((alpha)/2))*100),"%",sep=""))
ends2 <- c("SE","P.value")
CI <- c(lower=bvec[low],upper=bvec[up])
tr <- c(se, pv) 
res <- list(ci=CI,head=head, head2 =head2, ends2 =ends2, ends=ends, test.res = tr)
class(res)<-"bci"
res
}

print.bci <- function(x,digits= max(3, getOption("digits")), ...){
cat("\n")
cat(x$head,"\n")
rq<-structure(x$ci,names=x$ends)
print(rq,digits=digits)
cat("\n")
cat(x$head2,"\n")
rq2<-structure(x$test.res,names=x$ends2)
print(rq2,digits=digits)
cat("\n")
invisible(x)
}
