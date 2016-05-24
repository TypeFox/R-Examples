`predict.TAR` <-
function (object,n.ahead=1,n.sim=1000,...) 
{
res=NULL
p1=object$p1
p2=object$p2
d=object$d
xstart=rev(rev(object$y)[1:max(d,p1,p2)])
for (i in 1:n.sim) {
res=rbind(res, tar.sim(object,ntransient=0, n=n.ahead,xstart=xstart)$y)
}
colnames(res)=paste(1:n.ahead, "step",sep='-')
quan=apply(res,2,quantile, c(0.5,0.025,0.975))
invisible(list(pred.matrix=res, fit=quan[1,],pred.interval=quan[-1,]))
}

