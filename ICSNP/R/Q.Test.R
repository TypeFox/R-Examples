`Q.Test`<-function(X,score)
{
    n<-dim(X)[1]
    p<-dim(X)[2]
    T<-n^-0.5*colSums(sign(X)*switch(score,
              "sign"=1,
              "rank"=apply(abs(X),2,rank)/(n+1),
              "normal"=qnorm(0.5+apply(abs(X),2,rank)/(2*(n+1)))
              ))
    Q<-sum(T^2)*switch(score,"rank"=3,1)
    p.value<-1-pchisq(Q,p)

    list(test.statistic=Q,p.value=p.value)
}
