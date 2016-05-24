print.bict <-
function(x,...){

hrs<-round(x$time%/%3600,0)
mins<-round((x$time%%3600)%/%60,0)
secs<-round((x$time%%3600)%%60,0)
hrs<-ifelse(hrs<10,paste("0",hrs,sep=""),hrs)
mins<-ifelse(mins<10,paste("0",mins,sep=""),mins)
secs<-ifelse(secs<10,paste("0",secs,sep=""),secs)

priortypes<-c("UIP","SBH")

cat("Number of cells in table =",length(x$maximal.mod$y),"\n")
cat("\n")
cat("Maximal model =\n")
print(x$maximal.mod$formula)
cat("\n")
cat("Number of log-linear parameters in maximal model =",dim(x$maximal.mod$x)[2],"\n")
cat("\n")
cat("Number of MCMC iterations =",length(x$MODEL),"\n")
cat("\n")
cat("Computer time for MCMC =",paste(hrs,":",mins,":",secs,sep=""),"\n")
cat("\n")
cat("Prior distribution for log-linear parameters =",priortypes[x$priornum],"\n")
cat("\n")
cat("Number of missing cells =",length(x$missing1),"\n")
cat("\n")
cat("Number of censored cells =",length(x$missing2),"\n")

}
