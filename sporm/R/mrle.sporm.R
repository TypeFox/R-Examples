mrle.sporm <-
function(x, y, theta=1, p=rep(1/(length(x)+length(y)), length(x)+length(y)), 
    tol=1e-7, maxit=50){
    m<-length(x)
    n<-length(y)
    N<-m+n
    z<-c(y,x)
    r<-rank(z, ties.method = "max")[1:n]
    del<-1
    it<-0
    res<-grad.hessinv(theta,p, r)
    del<-sum(abs(res$H))+abs(sum(p)-1)
    while(del>tol && it<maxit){
        theta.p<-c(theta, p)-try(res$Ainv%*%res$H, TRUE)
        res<-grad.hessinv(theta.p[1], theta.p[-1], r)
        del<-sum(abs(res$H))+abs(theta-theta.p[1])+abs(sum(theta.p[-1])-1)
        theta<-theta.p[1]
        p<-theta.p[-1]
        #cat(it+1, "del=", del, theta, "sum of p=", sum(p),"\n")
        it<-it+1
    }
    if(del>tol || theta<=0 ||abs(sum(p)-1)>1e-2 || any(p<0) || any(p>1)) stop("Not Convergent\n")
    temp<-NULL
    temp1<-NULL
    for(j in 1:n){
        temp[j]<-sum(p[1:r[j]])
        temp1[j]<-ifelse(r[j]==1,0,sum(p[1:(r[j]-1)]))
    }
    ell<-n*log(theta)+sum(log(p))-sum(log(1+(theta-1)*temp))-sum(log(1+(theta-1)*temp1))
    list(theta=theta, p=p, ell=ell, it=it, del=del)
}
