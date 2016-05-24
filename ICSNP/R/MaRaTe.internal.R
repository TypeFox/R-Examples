`MaRaTe.internal`<-function(X,Y=NULL,test)
{
 n<-dim(X)[1]
 p<-dim(X)[2]
 Xsigns<-sign(X)
 Xranks<-apply(abs(X),2,rank)
 if(is.null(Y))  # one sample
 {
  scores<-switch(test,"sign"= Xsigns, "normal"=Xsigns*qnorm(0.5+Xranks/(2*(n+1))), "rank"=Xsigns*Xranks)
  A<-t(scores)%*%scores/switch(test,"rank"=(n+1)^2,1)
  S<-colSums(scores)/switch(test,"rank"=(n+1),1)
  test.statistic<-as.numeric(t(S)%*%solve(A)%*%S)
  p.value<-1-pchisq(test.statistic,p)
  return(list(test.statistic=test.statistic,p.value=p.value))
 }
 # else two samples
 m<-dim(Y)[1]
 
 XYdata<-rbind(X,Y)
 XYranks<-apply(XYdata,2,rank)
 E<-switch(test, "sign"={medians<-apply(XYdata,2,median); 
                          as.matrix(ifelse(XYdata<=medians,1,0));},
                  "rank"=XYranks/(n+m+1),
                  "normal"=qnorm(XYranks/(n+m+1)))
 E.X<-E[1:n,]
 E.Y<-E[(n+1):(n+m),]
 T.X<-mycolMeans((E.X))
 T.Y<-mycolMeans((E.Y))
 E.bar<-mycolMeans(E)
 W.x <- t(E.X)%*%E.X - n*(E.bar%*%t(E.bar))
 W.y <- t(E.Y)%*%E.Y - m*(E.bar%*%t(E.bar))
 W<-(1/(n+m))*(W.x+W.y)
 W.inv<-solve(W)
 L<-as.numeric(n*t(T.X-E.bar)%*%W.inv%*%(T.X-E.bar)+m*t(T.Y-E.bar)%*%W.inv%*%(T.Y-E.bar))
                        
 p.value<-1-pchisq(L,p)

 list(test.statistic=L,p.value=p.value,Cov=W,df=p)
}


mycolMeans<-function(X){
apply(as.data.frame(X),2,mean)}


MaRaTe.internal.csample<-function(X,g,scores="rank",...)
 {
 Xranks<-apply(X,2,rank)
 g.levels<-levels(g)
 
 scores<-match.arg(scores,c("sign","rank","normal"))
 N<-dim(X)[1]
 if (N!=length(g)) stop("number of rows of 'X' must correspond to the length of 'g'")
 p<-dim(X)[2]
 E<-switch(scores, "sign"={medians<-apply(X,2,median); 
                          as.matrix(ifelse(X<=medians,1,0));},
                  "rank"=Xranks/(N+1),
                  "normal"=qnorm(Xranks/(N+1)))
 E.g<-split(as.data.frame(E),g)
 
 n.g<-by(E[,1],g,length)
 if (min(n.g)<2 | is.na(min(n.g))) stop("for each factor level of 'g' at least two observations are needed")
 
 T.g<-by(E,g,mycolMeans)
 E.bar<-mycolMeans(E)
 W.g<-sapply(sapply(E.g,as.matrix,simplify=F),crossprod,simplify=F)
 W.bar<-sapply(n.g,"*",E.bar%*%t(E.bar),simplify=F)
 
 
 W.i <- matrix(unlist(lapply(W.g, t)), ncol=p, byrow=TRUE) - matrix(unlist(lapply(W.bar, t)), ncol=p, byrow=TRUE)
 W.i <- split(as.data.frame(W.i),rep(g.levels,each=p))
 
 W<- matrix(0,p,p)
 for (i in seq_along(g.levels))
    {
    W<-W+matrix(unlist(W.i[i]),p)
    }
 W<-W/N
 W.inv<-solve(W)
 
 L<-0
 diff.TE <- sapply(T.g,"-",E.bar,simplify=F) 
 
 for (i in seq_along(g.levels))
    {
    L<-L+n.g[i]* (matrix(unlist(diff.TE[i]),1,p) %*%W.inv%*% (matrix(unlist(diff.TE[i]),p)))
    }
 
 L<-as.numeric(L)
 dfs<-p*(length(g.levels)-1)                       
 p.value<-1-pchisq(L,dfs)

 list(test.statistic=L,p.value=p.value,Cov=W,df=dfs)
 }
