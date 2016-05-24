model.eval <- function (t, Xobs, Xmod, fileout, pdfout=F) {

# absolute error: individual, total and average
abs.err <- Xmod-Xobs
abs.err.tot <- sum (abs.err)
abs.err.avg <- abs.err.tot/length(Xobs)

# relative error in percent: individual, total and average
rel.err <- 100*(Xmod-Xobs)/Xobs
rel.err.tot <- sum(rel.err)
rel.err.avg <- rel.err.tot/length(Xobs)

# squared error: individual, total and average
sq.err <- (Xmod-Xobs)^2
sq.err.tot <- sum (sq.err)
sq.err.avg <- sq.err.tot/length(Xobs)

# rms error
rms.err <- sqrt(sq.err.avg)

# binding the results together and round to two decimal places
result <- round(cbind(t, Xobs, Xmod, abs.err, rel.err, sq.err)*100)/100 
result.tot <- round(c(abs.err.tot, rel.err.tot, sq.err.tot)*100)/100
result.avg <- round(c(abs.err.avg, rel.err.avg, sq.err.avg, rms.err)*100)/100
result.tot.avg <- c(result.tot, result.avg)

# labels for the results
label.result <- c("time", "Xobs", "Xmod", "abs.err", "rel.err", "sq.err")
label.result.tot.avg <- c("abs.tot", "rel.tot", "sq.tot", "abs.avg", "rel.avg",                          "sq.avg","rms")

# write file with results
fileouttxt <- paste(fileout,".txt",sep="")
write(t(label.result), fileouttxt, ncolumns=6)
write(format(t(result), nsmall=2), fileouttxt, ncolumns=6, append=T)
write(t(c(label.result.tot.avg)), fileouttxt, ncolumns=7, append=T)
write(format(t(c(result.tot,result.avg)),nsmall=2), fileouttxt, ncolumns=7, append=T)

# graph to visualize use four panels
if(pdfout== T) pdf(file=paste(fileout,".pdf",sep=""))
 mat <- matrix(1:4,2,2,byrow=T)
 layout(mat, widths=rep(7/2,4), heights=rep(7/2,4), TRUE)
 par(mar=c(4,4,1,.5),xaxs="r", xaxs="r")

# dynamics plot: model and observed vs time
plot(t,Xobs, ylim=c(0,max(Xobs)))
lines(t,Xmod)
legend (min(t), min(Xobs)+20,legend=c("Xobs", "Xmod"),
       pch=c(1,-1), lty=c(-1,1), merge=T) 

# graph the errors use horizontal lines for zero and averages
plot(t, abs.err); abline(h=0);abline(h=abs.err.avg, lty=c(2))
plot(t, rel.err); abline(h=0);abline(h=rel.err.avg, lty=c(2))
plot(t, sq.err); abline(h=0);abline(h=sq.err.avg, lty=c(2))

if(pdfout== F) {
 mat <- matrix(1:4,2,2,byrow=T)
 layout(mat, widths=rep(7/2,4), heights=rep(7/2,4), TRUE)
 par(mar=c(4,4,1,.5),xaxs="r", xaxs="r")
}

# scatter plot: observed vs model compared to unit slope line
plot(Xobs, Xmod, xlim=c(0,max(Xobs)))
lines(Xobs, Xobs)
legend (min(Xobs), max(Xmod),legend=c("Model", "Exact"),
       pch=c(1,-1), lty=c(-1,1), merge=T) 

if(pdfout== T) dev.off()

return(list(result=result, label.result.tot.avg=label.result.tot.avg, result.tot.avg=result.tot.avg))

}
