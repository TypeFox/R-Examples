# read-plot-seem-fileout-function.R 
# scan and import seem file_out output file to R

 read.plot.forsucc.out <- function(file.out.name, sens=F, pdfout=F){
 
# read the first three lines to get number of runs, of variables and data values 
 nvalues <- scan(file.out.name, what=c("",""), nlines=3)
 nruns <- as.numeric(nvalues[2])
 nvars=as.numeric(nvalues[6])
 ndata=as.numeric(nvalues[10]) 

# create array to store values of parameter
 val.par <- as.character(c(1:nruns))
 label.var <- scan(file.out.name, skip= 6, nlines=1, what=c(""), nmax=nvars+1)

label.var <- c("Years", "Sp1","Sp2","Sp3","Sp4", "Tot", "Sp1","Sp2","Sp3","Sp4", "Tot")
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
 for(i in c(1:4)){
 xvar <- x[,,i]
 # produce graph 
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab=label.var[i+1],lty=c(1:nruns), col=1)
  legend (1.05*max(t), max(xvar), legend=val.par, lty=c(1:nruns), col=1)
}

} else{
# produce graph

  plot.vars <- c(1:4)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Relative BA",lty=c(1:nplotvars), col=1)
  legend (1.05*max(t), max(xvar), legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1)

  plot.vars <- c(5)
  xvar <- x[,1,plot.vars] *(100/10000) # convert BA to m2 and sim area to ha
  nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="BA [m2/ha]",lty=c(1:nplotvars), col=1)
  legend (1.05*max(t), max(xvar), legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1)

  plot.vars <- c(6:9)
  xvar <- x[,1,plot.vars]
  nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Rel Density",lty=c(1:nplotvars), col=1)
  legend (1.05*max(t), max(xvar), legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1)

  plot.vars <- c(10)
  xvar <- x[,1,plot.vars] # convert to density ind/ha sim aea to ha but give as dens/100
  nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Density [stems/ha]/100",lty=c(1:nplotvars), col=1)
  legend (1.05*max(t), max(xvar), legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1)

  if(pdfout==F) win.graph()
}
  dev.off()
} else{

mat <- matrix(1:2,2,1,byrow=T)
nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")

if(sens==T){
 # select sens variable to plot
 for(i in c(1:4)){
 xvar <- x[,,i]
 # produce graph 
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab=label.var[i+1],lty=c(1:nruns), col=1)
  legend (1.05*max(t), max(xvar), legend=paste(lab.par,"=",val.par), lty=c(1:nruns), col=1)

   if(i%%2==0) {
   win.graph(); nf <- layout(mat, widths=rep(7,2), heights=rep(7/2,2), TRUE)
   par(mar=c(4,4,3,.5),xaxs="i",yaxs="i")
  }
 }
}

else{
# produce graph

  plot.vars <- c(1:4)
  xvar <- x[,1,plot.vars]; nplotvars <- length(plot.vars)
  matplot(t, xvar, type="l", ylim=c(0,max(xvar)), xlim=c(0,1.5*max(t)),
           xlab=label.var[1],ylab="Population & Nutrients",lty=c(1:nplotvars), col=1)
  legend (1.05*max(t), max(xvar), legend=label.var[plot.vars+1], lty=c(1:nplotvars), col=1)

}
}
 
 return (list(t=t, x=x, val.par=val.par, label.var=label.var))
} 

