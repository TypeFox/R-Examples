comp.graph <-
function(obj,alpha=0.05,xlab="time",ylab="X effect",all.terms=TRUE,variable) {
    se <- solve(obj$hess)%*%obj$h2%*%solve(obj$hess)
    dm <- ncol(obj$hess)
    sq <- seq(1,dm,by=ncol(obj$Ft))
    if(all.terms==TRUE){
        for(i in 1:ncol(obj$X)){
            f2 <- obj$theta[i,]%*%t(obj$Ft)
            Dy <- matrix(rep(0,ncol(obj$hess)*length(f2)),ncol=dm)
            Dy[1:nrow(obj$X),sq[i]:(sq[i]+ncol(obj$Ft)-1)] <- t(obj$theta[i,]*t(obj$Ft))
            vy2 <- Dy %*% se %*% t(Dy)
            B2se <-  sqrt(diag(abs(vy2)))
            ul <- f2+qnorm(1-alpha/2)*B2se
            ll <- f2-qnorm(1-alpha/2)*B2se
            plot(obj$time,f2,'l',xlab=xlab,ylab=ylab,lwd=3,ylim=c(min(ll),max(ul)))
            polygon(c(obj$time,obj$time[ncol(f2):1],obj$time[1]), c(ul, ll[ncol(f2):1], ul[1]), col = "gray",border=NA)
            lines(obj$time,f2,lwd=3,lty=2)}}
    if(all.terms==FALSE){
        dim <- variable
        f2 <- obj$theta[dim,]%*%t(obj$Ft)
        Dy <- matrix(rep(0,ncol(obj$hess)*length(f2)),ncol=ncol(obj$hess))
        Dy[1:nrow(obj$X),sq[dim]:(sq[dim]+ncol(obj$Ft)-1)] <- t(obj$theta[dim,]*t(obj$Ft))
        vy2 <- Dy %*% se %*% t(Dy)
        B2se <-  sqrt(diag(abs(vy2)))
        ul <- f2+qnorm(1-alpha/2)*B2se
        ll <- f2-qnorm(1-alpha/2)*B2se
        plot(obj$time,f2,'l',xlab=xlab,ylab=ylab,lwd=3,ylim=c(min(ll),max(ul)))
        polygon(c(obj$time,obj$time[ncol(f2):1],obj$time[1]), c(ul, ll[ncol(f2):1], ul[1]), col = "gray",border=NA)
        lines(obj$time,f2,lwd=3,lty=2)}
}
