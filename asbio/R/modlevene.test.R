modlevene.test<-function(x, groups){
medians<-sapply(split(x,groups),median,na.rm=TRUE)
resid.x<-abs(x-medians[groups])
res<-list()
res<-anova(lm(resid.x~groups))
class(res) <- "mltest"
res}

print.mltest <- function(x, ...){
df1 <- paste("df1 =", x$"Df"[1])
df2 <- paste("df2 =", x$"Df"[2])
F <-  paste("F =", round(x$"F value"[1], 5))
P <- paste("p-value =", round(x$"Pr(>F)"[1], 5))  
out <- paste(df1, df2, F, P, sep = ",  ")
cat("\n")
cat("Modified Levene's test of homogeneity of variances","\n")
cat("\n")
cat(out)
invisible(x)
}