# read-plot-seem-fileout-function.R 
# scan and import seem file_out output file to R

 read.plot.cerio.out <- function(file.out.name, sens=F, pdfout=F){
 
# read the first three lines to get number of runs, of variables and data values 
 nvalues <- scan(file.out.name, what=c("",""), nlines=3)
 nruns <- as.numeric(nvalues[2])
 nvars=as.numeric(nvalues[4])
 ndata=as.numeric(nvalues[6]) 

# create array to store values of parameter
 val.par <- as.character(c(1:nruns))
 label.var <- scan(file.out.name, skip= 6, nlines=1, what=c(""), nmax=nvars+1)

# extract the first column of one of the runs for the time stamps
 t <- matrix(scan(file.out.name,skip= 8,nlines=ndata), ncol=nvars+1, byrow=T)[,1]
# 3d array to store the results
 x <- structure(1:(nruns*ndata*nvars), dim=c(ndata,nruns,nvars))

# loop thru runs extract labels and variable values, bind to results
 for(i in 1:nruns) {
   valpart <- scan(file.out.name, skip= 7 + (i-1)*(ndata+1),nlines=1, 
                     what=c("","","","","",""))
   lab.par <- valpart[4]
   val.par[i] <- valpart[length(valpart)]

   x[,i,] <- matrix(scan(file.out.name, skip= 8 + (i-1)*(ndata+1),nlines=ndata ),
                      ncol=nvars+1, byrow=T)[,-1]
 }


if(pdfout==T){
pdf(paste(file.out.name,".pdf",sep=""))
mat <- matrix(1:2,2,1,byrow=T)
nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")
if(sens==T){
 # select sens variable to plot
 for(i in c(1:7)){
 xvar <- x[,,i]
 # produce graph 
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab=label.var[i+1],lty=c(1:nruns), col=1)
  legend ("topright", legend=paste(lab.par,"=",val.par), lty=c(1:nruns), col=1, cex=0.8)
}

} else{
# produce graph

  plot.vars <- c(1:4,7)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Population",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)

  plot.vars <- c(5:6)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Algae Concentration",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)

  if(pdfout==F) win.graph()
 
  plot.vars <- c(8:10)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Conditions",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)
}
  dev.off()
} else{

mat <- matrix(1:2,2,1,byrow=T)
nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")

if(sens==T){
 # select sens variable to plot
 for(i in c(1:7)){
 xvar <- x[,,i]
 # produce graph 
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab=label.var[i+1],lty=c(1:nruns), col=1)
  legend ("topright", legend=paste(lab.par,"=",val.par), lty=c(1:nruns), col=1, cex=0.8)

   if(i%%2==0) {
   win.graph(); nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
   par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")
  }
 }
}

else{
# produce graph

  plot.vars <- c(1:4,7)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Population",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)

  plot.vars <- c(5:6)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Algae Concentration",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)

   win.graph(); nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
   par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")

  plot.vars <- c(8:10)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Conditions",lty=c(1:nplotvars), col=1)
  legend ("topright", legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1, cex=0.8)
}
}
 
 return (list(t=t, x=x, val.par=val.par, label.var=label.var))
} 

