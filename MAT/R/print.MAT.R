print.MAT <-
function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\n")
stats.NIA<-c(summary(x$ni.administered),sd(x$ni.administered))
names(stats.NIA)<-c("MIN","Q1","Median","Mean","Q3","MAX","SD")
cat("  Number of items administered:\n")
print(stats.NIA,digits=2)
cat("\n")
cat("  Correlation between CAT and full bank thetas:\n")
Thetas<-data.frame(x$theta.Full,x$theta.CAT)
names(Thetas)<-c(paste("CAT",1:x$p),paste("Bank",1:x$p))
cor.Theta<-cor(Thetas)
print(cor.Theta,digits=3)
cat("\n")
}

