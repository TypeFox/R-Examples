`seg.control` <-
function(toll=.0001, it.max=10, display=FALSE, stop.if.error=TRUE, K=10, quant=FALSE, last=TRUE, maxit.glm=25, h=1,
    n.boot=20, size.boot=NULL, gap=FALSE, jt=FALSE, nonParam=TRUE, random=TRUE, powers=c(1,1), seed=NULL, fn.obj=NULL){
        list(toll=toll,it.max=it.max,visual=display,stop.if.error=stop.if.error,
            K=K,last=last,maxit.glm=maxit.glm,h=h,n.boot=n.boot, size.boot=size.boot, gap=gap, jt=jt,
            nonParam=nonParam, random=random, pow=powers, seed=seed, quant=quant, fn.obj=fn.obj)}

