ind.ctest<-function(X,index1,index2 = NULL,scores="rank",na.action=na.fail)
    {
    dname<-deparse(substitute(X))
    
    dindex1<-deparse(substitute(index1))
    
    if(is.null(index2)) {dindex2<-paste("-(",dindex1,")",sep = "")}
    else
    {dindex2<-paste("(",deparse(substitute(index2)),")",sep = "")}
    
    DNAME <- paste(dname,"[,(",dindex1,")]"," and ",dname,"[,",dindex2,"]",sep = "")
    
    if(!is.null(index2)) {if(any(index1 %in% index2)) 
    stop("'index1' and 'index2' cannot select the same columns")}
    
    
    p<-dim(X)[2]
    
    if (length(index1)==p) stop ("'index1' cannot select all or no columns of 'X'")
    
    X1<-as.matrix(X[,index1])
    p.ind1 <- dim(X1)[2]
    
    if(is.null(index2)) {X2<-as.matrix(X[,-index1])}
    else
    {X2<-as.matrix(X[,index2])}
    p.ind2<-dim(X2)[2]
    X<-cbind(X1,X2)
    X<-na.action(X)
    n<-dim(X)[1]
    
    if(!all(sapply(X, is.numeric))) stop("all selected columns must be numeric")
    X<-as.matrix(X)
    
    X.ranks <- apply(X,2,rank)
    scores <- match.arg(scores,c("sign","rank","normal"))
    
    scores1 <- function(x){  X.med<-apply(x,2,median)
                               Y<-sweep(x,2,X.med,"-")
                               apply(Y,2,sign)}
    scores2 <- function(x,n){(12/(n^2+1))^0.5*(x-(n+1)/2)}  
    
    E<-switch(scores, "sign"=scores1(X),
                  "rank"=apply(X.ranks,2,scores2,n),
                  "normal"=apply(X.ranks/(n+1),2,qnorm))  
    
    T<- t(E)%*%E/n
    
    W<-det(T)/(det(as.matrix(T[1:p.ind1,1:p.ind1]))*det(as.matrix(T[-(1:p.ind1),-(1:p.ind1)])))
    
    W.stat<- -n*log(W)
    p.value<- 1-pchisq(W.stat,p.ind1*p.ind2)
    
    STATISTIC<-W.stat
    names(STATISTIC)<-"W"
    
    Scores<-switch(scores, "sign"="Signs",
                  "rank"="Ranks",
                  "normal"="Normal Scores")  
    
    METHOD<-paste("Test of Independence Based on Marginal", Scores)
    PARAMETER<-p.ind1*p.ind2
    names(PARAMETER)<-c("df")
    PVAL<-p.value
    res<-list(method=METHOD,statistic=STATISTIC,data.name=DNAME,parameter=PARAMETER,p.value=PVAL)
    class(res)<-"htest"
    return(res)
    
    }
