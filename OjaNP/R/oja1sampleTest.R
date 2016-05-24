`oja1sampleTest` <-
function(X, mu=NULL, scores="sign", p = 1, method = "approximation", n.simu = 1000, na.action=na.fail, ...)
    {
    DNAME<-deparse(substitute(X))
    
    n<-dim(X)[1]
    k<-dim(X)[2]

    if (is.null(mu)) mu<-rep(0,k)
    else if (length(mu)!=k) stop("length of 'mu' must equal the number of columns of 'X'")
    
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)

    X<-sweep(X,2,mu)
    score<-match.arg(scores,c("sign","rank"))
    method <- match.arg(method,c("approximation","permutation"))
    
    CENTER<- rep(0,k)
    
    res1<-switch(score,
        "sign"={
               SIGNS <- ojaSign(X, x= NULL, center=CENTER, p=p, ...)
               SUM.SIGNS <- colSums(SIGNS)
               A.inv <- solve(crossprod(SIGNS)) 
               STATISTIC<- as.numeric((SUM.SIGNS) %*% A.inv %*% (SUM.SIGNS))
               names(STATISTIC)<-"Q.S"
               
               switch(method,
                             "approximation"={
                             PVAL <- 1-pchisq(as.numeric(STATISTIC),k)
                             PARAMETER<-k
                             names(PARAMETER)<-"df" 
                             }
                             ,
                             "permutation"={
                             rep.func.sign <- function(index,signs,A.inv)
                                {
                                SUM.SIGNS.2 <- colSums(index*signs)
                                as.vector((SUM.SIGNS.2) %*% A.inv %*% (SUM.SIGNS.2))
                                }
                             Q.simu <- replicate(n.simu, rep.func.sign(sample(c(1,-1), n, replace = TRUE), SIGNS, A.inv))
                             PVAL <- mean(Q.simu>STATISTIC)
                             PARAMETER<-n.simu
                             names(PARAMETER)<-"replications"
                             }
                             )
               METHOD<-"OJA 1 SAMPLE SIGN TEST"
               RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD)
               
               RVAL}
        ,
        "rank"={
               SIGNEDRANKS <- ojaSignedRank(X, x=NULL, p=p, ...)
               SUM.SR <- colSums(SIGNEDRANKS)
               B.inv <- solve(crossprod(SIGNEDRANKS)) 
               STATISTIC<- as.numeric((SUM.SR) %*% B.inv %*% (SUM.SR))
               names(STATISTIC)<-"Q.R"
              
               switch(method,
                             "approximation"={
                             PVAL <- 1-pchisq(as.numeric(STATISTIC),k)
                             PARAMETER<-k
                             names(PARAMETER)<-"df" 
                             }
                             ,
                             "permutation"={
                             rep.func.rank <- function(index,signs,B.inv)
                                {
                                SUM.SR.2 <- colSums(index*signs)
                                as.vector((SUM.SR.2) %*% B.inv %*% (SUM.SR.2))
                                }
                             Q.simu <- replicate(n.simu, rep.func.rank(sample(c(1,-1), n, replace = TRUE), SIGNEDRANKS, B.inv))
                             PVAL <- mean(Q.simu>STATISTIC)
                             PARAMETER<-n.simu
                             names(PARAMETER)<-"replications"
                             }
                             )
               METHOD<-"OJA 1 SAMPLE SIGNED RANK TEST"
               RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD)
               
               RVAL
               }
        )
    
    ALTERNATIVE="two.sided"
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    names(NVAL)<-"location"
    res<-c(res1,list(data.name=DNAME,parameter=PARAMETER,alternative=ALTERNATIVE,null.value=NVAL))
    class(res)<-"htest"
    return(res)
    }
