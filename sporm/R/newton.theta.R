newton.theta <-
function(y, theta0=1, maxit = 100, eps = 1e-10){
    n<-length(y)
    score<-function(theta) n/theta-2*sum(y/(1+(theta-1)*y))
    dscore<-function(theta) -n/theta^2+2*sum(y^2/(1+(theta-1)*y)^2)
    del<-abs(score(theta0))
    it<-1
    while(del>eps & it<maxit){
        theta<-theta0-score(theta0)/dscore(theta0)
        del<-abs(score(theta0))+abs(theta0-theta)
        it<-it+1
        theta0<-theta
    }
    theta
}

