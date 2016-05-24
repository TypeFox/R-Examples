
mycolMeans<-function(X){
apply(as.data.frame(X),2,mean)}



ojaCsampleTestFormula<-function(X,g,scores="rank",method="approximation", n.simu=1000, p = 1, center = "ojaMedian", ...)
 {
 g.levels<-levels(g)
 
 scores<-match.arg(scores,c("sign","rank"))
 N<-dim(X)[1]
 if (N!=length(g)) stop("number of rows of 'X' must correspond to the length of 'g'")
 k<-dim(X)[2]
 
 
 
 E<-switch(scores, "sign"= ojaSign(X, center=center, p=p, ...),
                   "rank"= ojaRank(X, p=p, ...))
 E.g<-split(as.data.frame(E),g)
 
 n.g<-by(E[,1],g,length)
 if (min(n.g)<2 | is.na(min(n.g))) stop("for each factor level of 'g' at least two observations are needed")
 
 T.g<-by(E,g,mycolMeans)
 Bmat <- crossprod(sweep(t(sapply(T.g, cbind)),1,sqrt(n.g),"*"))
 Dmat <- crossprod(E)/N
 
 L <- sum(diag(solve(Dmat,Bmat)))

 STATISTIC<-as.numeric(L)
 
 switch(method,
            "approximation"={
            PARAMETER<-k*(length(g.levels)-1) 
            PVAL <- 1-pchisq(STATISTIC, PARAMETER)
            names(PARAMETER)<-"df" 
            }
            ,
            "permutation"={
            rep.func.sign <- function(index,g,score,Dmat.inv,n.g)
                                {
                                g2 <- g[index]
                                T.g2<-by(score,g2,mycolMeans)
                                Bmat2 <- crossprod(sweep(t(sapply(T.g2, cbind)),1,sqrt(n.g),"*"))
                                L2 <- sum(diag(crossprod(Dmat.inv,Bmat2)))
                                as.numeric(L2)
                                }
            Dmat.inv <- solve(Dmat)
            Q.simu <- replicate(n.simu, rep.func.sign(sample(1:N), g, E,  Dmat.inv,n.g))
            PVAL <- mean(Q.simu>STATISTIC)
            PARAMETER<-n.simu
            names(PARAMETER)<-"replications"
                             }
                             )
 names(STATISTIC)<-ifelse(scores=="sign", "Q.S", "Q.R")
 list(statistic=STATISTIC,p.value=PVAL,df=PARAMETER)
 }
