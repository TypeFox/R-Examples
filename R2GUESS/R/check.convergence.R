check.convergence <-
function(x,which = c(1L:2L),nsplit=10,nbloc=4,ask = prod(par("mfcol")) < length(which) && dev.interactive()){

show <- rep(FALSE, 2)
show[which] <- TRUE

if(is.null(x$label.Y)) Pheno <- paste("Y",1:x$q,sep="",collapse="_") else Pheno <- paste("Y: ",paste(x$label.Y,collapse="/"),sep="") 

NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl.txt",sep="_"))
NameGHistory <- file.path(x$path.output, paste(x$root.file.output,"output_g_history.txt",sep="_"))
NameModelSize <- file.path(x$path.out, paste(x$root.file.output,"output_model_size_history.txt",sep="_"))
NameCondPost <- file.path(x$path.out, paste(x$root.file.output,"output_log_cond_post_prob_history.txt",sep="_"))
NameTemp <- file.path(x$path.out, paste(x$root.file.output,"output_temperature_history.txt",sep="_"))
Marg <- read.table(NameMarg,header=TRUE)
NameModelSize <- file.path(x$path.out, paste(x$root.file.output,"output_model_size_history.txt",sep="_"))
modelSize <- read.table(NameModelSize,header=TRUE)
modelSize <- modelSize[,2]
NameConv <- file.path(x$path.out, paste(x$root.file.output,"output_models_history.txt",sep="_"))

var.ord.MPI <- order(Marg[,2],decreasing=TRUE) ##representation of the variable with high MPPI


temp <- scan(NameConv,skip=2)
#####identify data without covariates
seq3 <- NULL
seq1 <- seq(1,2*x$nsweep,by=2)
seq2 <- seq(2,2*x$nsweep,by=2)
seq3[seq1] <- 4
seq3[seq2] <- modelSize
indix <- rep(rep(c(T,F),x$nsweep),times=seq3)
###########################################

my.data <- matrix(temp[indix],ncol=4,nrow=x$nsweep,byrow=TRUE)
colnames(my.data) <- c("nswwep","modelsize","log_marg","log_cond_post")

######################################################################################################
### log marginal
#####################################################################################################
mid <- floor((x$nsweep-x$burn.in)/2)

dens.all <- density(my.data[-c(1:x$burn.in),"log_marg"])     #   ,main="log-marginal",col="red")
dens.first <- density(my.data[(x$burn.in:(x$burn.in+mid)),"log_marg"])#  ,main="",col="blue",add=TRUE)
dens.last <- density(my.data[(x$burn.in+mid+1):x$nsweep,"log_marg"])#  ,main="",col="black",add=TRUE)
ymin <- min(dens.all$y,dens.first$y,dens.last$y)
ymax <- max(dens.all$y,dens.first$y,dens.last$y)
xmin <- min(dens.all$x,dens.first$x,dens.last$x)
xmax <- max(dens.all$x,dens.first$x,dens.last$x)

if(FALSE){
title.1 <- expression(paste("Log marginal likelihood distribution: log ",P(Y~group("|",list(gamma,tau),""))))
plot(dens.all,main=title.1,col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
par(new=TRUE)
plot(dens.first,main="",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
par(new=TRUE)
plot(dens.last,main="",col="green",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
text.legend <-c(paste("ALL= [",x$burn.in+1,":",x$nsweep,"]",sep=""),paste("First half= [",x$burn.in+1,":",(x$burn.in+mid),"]",sep=""),paste("Last half= [",x$burn.in+mid+1,":",x$nsweep,"]",sep=""))
legend("topleft",title="sweep",legend=text.legend,col=1:3,lty=1,text.col=1:3, bg = 'gray90')
}


######################################################################################################
### log cond post
#####################################################################################################


dens.all <- density(my.data[-c(1:x$burn.in),"log_cond_post"])     #   ,main="log-cond-post",col="red")
dens.first <- density(my.data[(x$burn.in:(x$burn.in+mid)),"log_cond_post"])#  ,main="",col="blue",add=TRUE)
dens.last <- density(my.data[(x$burn.in+mid+1):x$nsweep,"log_cond_post"])#  ,main="",col="black",add=TRUE)

ymin <- min(dens.all$y,dens.first$y,dens.last$y)
ymax <- max(dens.all$y,dens.first$y,dens.last$y)
xmin <- min(dens.all$x,dens.first$x,dens.last$x)
xmax <- max(dens.all$x,dens.first$x,dens.last$x)

if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }


title.2 <- expression(paste("Log Posterior Distribution: log ",P(list(gamma,tau)~group("|",Y,""))))

if (show[1L]) {
dev.hold()
plot(dens.all,main=title.2,col="black",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
par(new=TRUE)
plot(dens.first,main="",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
par(new=TRUE)
plot(dens.last,main="",col="green",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
text.legend <-c(paste("ALL= [",x$burn.in+1,":",x$nsweep,"]",sep=""),paste("First half=[",x$burn.in+1,":",(x$burn.in+mid),"]",sep=""),paste("Last half= [",x$burn.in+mid+1,":",x$nsweep,"]",sep=""))
legend("topleft",title="sweep",legend=text.legend,col=1:3,lty=1,text.col=1:3, bg = 'gray90')

dev.flush()
}
#################################################################################################

###nsplit number of split of the sweep
mid <- floor((x$nsweep-x$burn.in)/nsplit)
text.legend <- NULL
ymax <- NULL
xmin <- NULL
xmax <-NULL
list.dens <- NULL
for (i in 1:nbloc){
indice <- (x$nsweep:(x$nsweep-i*mid+1))
dens <- density(my.data[indice,"log_cond_post"])
ymax <- max(ymax,dens$y)
xmin <- min(xmin,dens$x)
xmax <- max(xmax,dens$x)
list.dens <- c(list.dens,list(dens))
}

if (show[2L]) {
dev.hold()
for (i in 1:nbloc){
indice <- (x$nsweep:(x$nsweep-i*mid+1))
plot(list.dens[[i]],main=title.2,col=i,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="")
text.legend <-c(text.legend,paste("set ",i," = [",min(indice),":",max(indice),"]",sep=""))
par(new=TRUE)  
}
legend("topleft",title="moving window",legend=text.legend,col=1:nbloc,lty=1,text.col=1:nbloc, bg = 'gray90')
dev.flush()
}
#################################################################################################
}



