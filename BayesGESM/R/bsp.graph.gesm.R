bsp.graph.gesm <-
function(object,which,var, xlab,ylab,main){

reem <- function(aa,b){
                ag <- aa
                ag <- sub("(", "", ag,fixed=TRUE)
                ag <- sub(")", "", ag,fixed=TRUE)
                ag <- sub(b, "", ag,fixed=TRUE)
                ag <- strsplit(ag, ",")
                ag <- ag[[1]][1]
                ag
        }

if(missingArg(which))
stop("The value of the which argument is missing!!",call.=FALSE)

if(which!=1 & which!=2)
stop("value of the which argument is invalid!!",call.=FALSE)

chains <- object$chains     
p <- object$p
q <- object$q
ks <- object$ks
ks2 <- object$ks2
nu0 <- object$nu0

quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

n <- length(object$y)
R <- nrow(chains)

if(which==1){
if(sum(ks) == 0) stop("A nonparametric function was not specified in the location model", call.=FALSE)

nps <- colnames(object$nps)
nvar <- grepl(var,nps,fixed=TRUE)
if(sum(nvar) == 0) stop("There is not the nonparametric effect requested by the user!!",call.=FALSE)
nvar <- min(seq(1,length(nps),by=1)[nvar==1])

if(missingArg(xlab) || !is.character(xlab))xlab <- reem(nps[nvar],"bsp")
if(missingArg(ylab) || !is.character(ylab))ylab <- nps[nvar]
if(missingArg(main) || !is.character(main))main <- " "

lims <-  c(0,cumsum(ks))
cov <- object$nps[,nvar]
B <- object$B[,(lims[nvar]+1):(lims[nvar+1])]
fs <- matrix(0,n,R)
for(i in 1:R) fs[,i] <- B%*%chains[i,(lims[nvar]+1+p):(lims[nvar+1]+p)]  
li <- apply(fs, 1, quant005)
m <- apply(fs, 1, mean)
ls <- apply(fs, 1, quant095)
id <- sort(cov,index=TRUE)$ix
plot(cov[id], li[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], m[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=1, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], ls[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
}

if(which==2){
if(sum(ks2) == 0) stop("A nonparametric function was not specified in the location model", call.=FALSE)
if(sum(ks)==0) length(ks) <- 0
nps2 <- colnames(object$nps2)
nvar <- grepl(var,nps2,fixed=TRUE)
if(sum(nvar) == 0) stop("There is not the nonparametric effect requested by the user!!",call.=FALSE)
nvar <- min(seq(1,length(nps2),by=1)[nvar==1])

if(missingArg(xlab) || !is.character(xlab))xlab <- reem(nps2[nvar],"bsp")
if(missingArg(ylab) || !is.character(ylab))ylab <- nps2[nvar]
if(missingArg(main) || !is.character(main))main <- " "

lims <-  c(0,cumsum(ks2))
cov <- object$nps2[,nvar]
D <- object$D[,(lims[nvar]+1):(lims[nvar+1])]
fs <- matrix(0,n,R)
for(i in 1:R) fs[,i] <- D%*%chains[i,(lims[nvar]+1+p+q+sum(ks)+length(ks)):(lims[nvar+1]+p+q+sum(ks)+length(ks))]  
li <- apply(fs, 1, quant005)
m <- apply(fs, 1, mean)
ls <- apply(fs, 1, quant095)
id <- sort(cov,index=TRUE)$ix
plot(cov[id], li[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], m[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=1, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], ls[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
}

}
