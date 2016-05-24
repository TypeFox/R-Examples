print.ridit <-
function (x,g,...) {
cat("\n")
cat("Ridit Analysis:")
cat("\n")
cat("\n")
#cat("Group:", "\t", paste(names(x$MeanRidit),"\t",sep="") ,"\n",sep="")
#cat("Ridit:", "\t", paste(round(x$MeanRidit,3),"\t",sep="") ,"\n",sep="")
### replacement for above two lines
m=max(nchar(names(x$MeanRidit)))
cutpoint=40
if(m>cutpoint) m=cutpoint
cat("Group","\t",format("Label",width=m), "\t", "Mean Ridit","\n",sep="")
cat("-----","\t",format("-----",width=m), "\t", "----------","\n",sep="")
for (k in 1:length(x$MeanRidit))
cat(k,"\t",format(substr(names(x$MeanRidit)[k],start=1,stop=m),width=m),"\t",round(x$MeanRidit,4)[k] ,"\n",sep="")
###
cat("\n")
cat(x$msg,"\n")
cat("chi-squared = ",round(x$TestStatistic,4),sep="")
cat(", df = ",x$df,sep="")
cat(", p-value = ",format(x$Sig,digits=4),sep="")
cat("\n")
cat("\n")
}
