T.vdw.loc <- function(mu,y)
    {
    n<-length(y)
    z<- y-mu
    T.mu <- -sum(sign(z)* qnorm(rank(abs(z))/(2*(n+1))+0.5))
    return(T.mu)
    }

vdw.loc <- function(x, int.diff=10, maxiter=1000, na.action=na.fail)
    {
     if (!is.vector(x)) stop("'x' must be a numeric vector")
    x<-na.action(x)
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    sums<-pair.sum(as.matrix(x))
    y.sorted<-sort(c(x,sums/2))
    m<-length(y.sorted)
    iter<-0
    int.lower<-1
    int.upper<-m
    res<-NULL
    while(iter<maxiter)
        {
        if (int.upper-int.lower<=int.diff) {
                           crit<- apply(matrix(y.sorted[(int.lower):int.upper]),1,T.vdw.loc,y=x)
                           #print(cbind(y.sorted[(int.lower):int.upper],crit)) 
                           index.min<-which.min(abs(crit))
                           crit.min <- crit[index.min]
                           index.min2 <- ifelse(crit.min<0, index.min+1, index.min-1) 
                           crit.min2 <- crit[index.min2]
                           x1 <-y.sorted[index.min+int.lower-1]
                           x2 <-y.sorted[index.min2+int.lower-1]
                           if(x1==x2) {res <- x1} else{
                           Mid<-data.frame(x=c(x1,x2),crit=c(crit.min,crit.min2))
                           perc<-with(Mid, abs(crit[1])/ (max(crit)-min(crit)))
                           res<-with(Mid, x[1] + sign(crit[2]) * perc* (max(x)-min(x)))}
                           break} else
        {int.middle<- floor((int.lower +int.upper)/2)
        Mid<-T.vdw.loc(y.sorted[int.middle], x)
        if (Mid <0) {int.lower<-int.middle} else {int.upper<-int.middle}}

        iter <- iter+1
        }
    #if (is.null(res)) stop("max.iter reached")
    #return(list(est=res,criterion=Mid))
    return(res)
    }
