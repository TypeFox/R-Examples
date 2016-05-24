icc <-
function(data)
{
score <- as.matrix(na.omit(data))
n <- dim(score)[1]
p <- dim(score)[2]
data2 <- matrix(ncol=3,nrow=p*n)
attr(score,"dim") <- c(p*n,1)
data2[,1] <- score
subject <- as.factor(rep(1:n,p))
rater <- as.factor(rep(1:p,each=n))
data2[,2] <- subject
data2[,3] <- rater
ms <- anova(lm(score~subject+rater))[[3]]
names(ms) <- NULL
v.s <- (ms[1]-ms[3])/p
v.r <- (ms[2]-ms[3])/n
res <- ms[3]
icc.a <- v.s/(v.s+v.r+res)
icc.c <- v.s/(v.s+res)
result <- list("nb.subjects"=n,"nb.raters"=p,"subject.variance"=v.s,"rater.variance"=v.r,"residual"=res,"icc.consistency"=icc.c,"icc.agreement"=icc.a)
result
}
