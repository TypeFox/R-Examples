plot.ESS <- function(x,...){
  args <- list(...)
  range.SF <- args["range.SF"]
  range.MS <- args["range.MS"]
  range.LP <- args["range.LP"]

  if(is.null(x$label.Y)) Pheno <- paste("Y",1:x$q,sep="",collapse="_") else Pheno <- paste("Y: ",paste(x$label.Y,collapse="/"),sep="") 

NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl.txt",sep="_"))
NameGHistory <- file.path(x$path.output, paste(x$root.file.output,"output_g_history.txt",sep="_"))
NameModelSize <- file.path(x$path.out, paste(x$root.file.output,"output_model_size_history.txt",sep="_"))
NameCondPost <- file.path(x$path.out, paste(x$root.file.output,"output_log_cond_post_prob_history.txt",sep="_"))
NameTemp <- file.path(x$path.out, paste(x$root.file.output,"output_temperature_history.txt",sep="_"))
Marg <- read.table(NameMarg,header=TRUE)



LogCondPost <-read.table(NameCondPost,header=TRUE)
gHistory <- read.table(NameGHistory,header=TRUE)
ModelSize <- read.table(NameModelSize,header=TRUE)
Temperature <- read.table(NameTemp,skip=1)

colnames(Temperature) <- c("Sweep",paste("Chain",1:x$nb.chain,sep="_"))
#jpeg(OutnamePlots,height=1000,width=1500,units='px',quality=80)

####Plot Dashboard:
#BurnIn <- max(10000,min(gHistory[,1]))
BurnIn <- x$burn.in
if(x$nb.chain>5)legCol <- colors()[sample(1:600,x$nb.chain)]else{
legCol <- c('blue','red','green','black')
legCol <- legCol[1:x$nb.chain]}


legTxt <- paste("Chain",1:x$nb.chain,sep=" ")
par(mfrow=c(2,2))
if(!is.null(range.SF[[1]])){
  #print(range.SF[1])
  ymin <- range.SF[[1]][1]
  ymax <- range.SF[[1]][2]
  plot(gHistory[-c(1:BurnIn),1],type='l',gHistory[-c(1:BurnIn),2],xlab='Sweep',ylab='Shrinkage factor',ylim=c(ymin,ymax),col='blue',lty=1,main=Pheno,cex.main=1)
  }else{
    plot(gHistory[-c(1:BurnIn),1],type='l',gHistory[-c(1:BurnIn),2],xlab='Sweep',ylab='Shrinkage factor',col='blue',lty=1,main=Pheno,cex.main=1)
    
  }

xmin <- min(ModelSize$Sweep)
xmax <- max(ModelSize$Sweep)
temp <- paste("min(ModelSize$Chain_",1:x$nb.chain,")",sep="")
temp1 <-paste(temp,collapse=",")
temp2 <-paste("min(",temp1,")",sep="")
ymin <- eval(parse(text=temp2))
temp <- paste("max(ModelSize$Chain_",1:x$nb.chain,")",sep="")
temp1 <-paste(temp,collapse=",")
temp2 <-paste("max(",temp1,")",sep="")
ymax <- eval(parse(text=temp2))

if(!is.null(range.MS[[1]])){
  ymin <- range.MS[[1]][1]
  ymax <- range.MS[[1]][2]
}
plot(ModelSize$Sweep,type='l',eval(parse(text=paste("ModelSize$Chain_",x$nb.chain,sep=""))),xlab='Sweep',ylab='Model Size',col=legCol[x$nb.chain],lty=1,xlim=c(xmin,xmax),ylim=c(ymin,ymax+1))

for (i in (x$nb.chain-1):1){
  lines(ModelSize$Sweep,eval(parse(text=paste("ModelSize$Chain_",i,sep=""))),col=legCol[i],lty=1)
}

xleg <- xmin+(xmax-xmin)*0.1
yleg <- ymin+(ymax-ymin)*0.8

if(!is.null(range.LP[[1]])){
  ymin <- range.LP[[1]][1]
  ymax <- range.LP[[1]][2]
}else{
ymin <- min(LogCondPost[,-1])
ymax <- max(LogCondPost[,-1])
}
     
plot(LogCondPost$Sweep,type='l',eval(parse(text=paste("LogCondPost$Chain_",x$nb.chain,sep=""))),xlab='Sweep',ylab='Log posterior',col=legCol[x$nb.chain],lty=1,main="Target distribution",ylim=c(ymin,ymax))
for (i in (x$nb.chain-1):1){
lines(LogCondPost$Sweep,eval(parse(text=paste("LogCondPost$Chain_",i,sep=""))),col=legCol[i],lty=1)
}

#Temp
xmin <- min(Temperature$Sweep)
xmax <- max(Temperature$Sweep)

temp <- paste("min(Temperature$Chain_",1:x$nb.chain,")",sep="")
temp1 <-paste(temp,collapse=",")
temp2 <-paste("min(",temp1,")",sep="")
ymin <- eval(parse(text=temp2))
temp <- paste("max(Temperature$Chain_",1:x$nb.chain,")",sep="")
temp1 <-paste(temp,collapse=",")
temp2 <-paste("max(",temp1,")",sep="")
ymax <- eval(parse(text=temp2))

plot(Temperature$Sweep,type='l',Temperature$Chain_1,xlab='Sweep',ylab='Temperature',col=legCol[1],lty=1,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
for (i in 2:x$nb.chain){
lines(Temperature$Sweep,eval(parse(text=paste("Temperature$Chain_",i,sep=""))),col=legCol[i],lty=1)
}
xleg <- xmin+(xmax-xmin)*0.5
yleg <- ymin+(ymax-ymin)*0.8
legend(xleg,yleg,col=legCol,legend=legTxt,lty=1,lwd=1,cex=1)
#dev.off()
}
######################################################################################################
