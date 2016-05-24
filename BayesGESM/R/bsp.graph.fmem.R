bsp.graph.fmem <-
function(object,which,xlab,ylab,main){

reem <- function(aa,b){
                ag <- aa
                ag <- sub("(", "", ag,fixed=TRUE)
                ag <- sub(")", "", ag,fixed=TRUE)
                ag <- sub(b, "", ag,fixed=TRUE)
                ag <- strsplit(ag, ",")
                ag <- ag[[1]][1]
                ag
        }

ks <- object$ks
homo <- object$homo
if(which<1 || which>length(ks) || which!=floor(which)) stop(paste("The which argument must be an integer value between in 1 and ",length(ks)))
if(sum(ks)==0) stop("A nonparametric function was not specified in the model")

nvar <- reem(colnames(object$nps)[which],"bsp")

if(missingArg(xlab) || !is.character(xlab))xlab <- nvar
if(missingArg(ylab) || !is.character(ylab))ylab <- paste("f(",nvar,")")
if(missingArg(main) || !is.character(main))main <- " "

quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

n <- length(object$y)
chains <- object$chains
R <- nrow(chains)

cov <- object$nps[,which]
lims <-  c(0,cumsum(ks))
B <- object$B[,(lims[which]+1):(lims[which+1])]
p <- object$p
q <- object$q
fs <- matrix(0,n,R)
for(i in 1:R){
if(homo==1) fs[,i] <- B%*%chains[i,(p+3*q+2+lims[which]):(p+3*q+2+lims[which]+ks[which]-1)]
else fs[,i] <- B%*%chains[i,(p+3*q+1+lims[which]):(p+3*q+lims[which]+ks[which])]
}
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
