Ell.Theta <-
function(x, y, theta.hat, p.hat, theta1, theta2, n.theta=40, tol=1e-7, maxit=500){
    m<-length(x)
    n<-length(y)
    N<-m+n
    z<-c(y,x)
    r<-rank(z, ties.method = "max")[1:n]
    z<-sort(z)
    Theta1<-seq(theta1, theta.hat, len=n.theta)
    Theta2<-seq(theta.hat, theta2, len=n.theta)
    ell<-NULL
    p<-p.hat
    theta<-theta.hat
    temp<-NULL
    temp1<-NULL
    for(j in 1:n){
        temp[j]<-sum(p[1:r[j]])
        temp1[j]<-ifelse(r[j]==1,0,sum(p[1:(r[j]-1)]))
    }
    ell[n.theta]<-n*log(theta)+sum(log(p))-sum(log(1+(theta-1)*temp))-sum(log(1+(theta-1)*temp1))
    for( i in n.theta:2){
        theta<-Theta1[i-1]
        res<-ell.theta(theta, p, r, tol, maxit)
        ell[i-1]<-res$ell
        p<-res$p
    }
    p<-p.hat
    for( i in 2:n.theta){
        theta<-Theta2[i]
        res<-ell.theta(theta, p, r, tol, maxit)
        ell[n.theta+i-1]<-res$ell
        p<-res$p
    }
    list(ell=ell, Theta=c(Theta1, Theta2[-1]))
}

