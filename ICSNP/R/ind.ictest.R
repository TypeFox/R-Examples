ind.ictest<-function(X,index1,index2 = NULL,scores="rank",method = "approximation", n.simu = 1000,...,na.action=na.fail)
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
    
    scores <- match.arg(scores,c("sign","rank","normal"))
    method <- match.arg(method,c("approximation","simulation","permutation"))
    
    X1<-as.matrix(X[,index1])
    p.ind1 <- dim(X1)[2]
    
    if(is.null(index2)) {X2<-as.matrix(X[,-index1])}
    else
    {X2<-as.matrix(X[,index2])}
    p.ind2<-dim(X2)[2]
    
    X<-data.frame(X1=I(X1),X2=I(X2))
    X<-na.action(X)
    n<-dim(X)[1]
    
    if(!all(sapply(X, is.numeric))) stop("all selected columns must be numeric")
    
    
    if (p.ind1==1)
    {ic.X1 <- X$X1}
    else
    {ic.X1 <- ics.components(ics(X$X1,...))}
    
    if (p.ind2==1)
    {ic.X2 <- X$X2}
    else
    {ic.X2 <- ics.components(ics(X$X2,...))}
    
    switch(scores, {
                   "sign"=
                   loc.ic.X1 <- apply(ic.X1,2,median)
                   loc.ic.X2 <- apply(ic.X2,2,median)
                   },
                  "rank"= {
                  loc.ic.X1 <- apply(ic.X1,2,hl.loc)
                  loc.ic.X2 <- apply(ic.X2,2,hl.loc)
                  },
                  "normal"= {
                  loc.ic.X1 <- apply(ic.X1,2,vdw.loc)
                  loc.ic.X2 <- apply(ic.X2,2,vdw.loc)
                  })  
    
    
    
    Z.X1 <- sweep(ic.X1,2,loc.ic.X1,"-")
    Z.X2 <- sweep(ic.X2,2,loc.ic.X2,"-")
    
    #Zs<-cbind(Z.X1,Z.X2)
    
    Z1.signs <- apply(Z.X1,2,sign)
    Z1.ranks <- apply(abs(Z.X1),2,rank)
    
    Z2.signs <- apply(Z.X2,2,sign)
    Z2.ranks <- apply(abs(Z.X2),2,rank)
    
    C.stat<-switch(scores, "sign"=t(Z1.signs) %*% Z2.signs/n,
                  "rank"= {t(Z1.signs*Z1.ranks)%*%(Z2.signs*Z2.ranks) * 3/(n*(n+1)^2)},
                  "normal"= {
                  t(Z1.signs*apply((Z1.ranks/(n+1)+1)/2,2,qnorm))%*%(Z2.signs*apply((Z2.ranks/(n+1)+1)/2,2,qnorm))/n})  
    
    
    
    Q.stat<- n*frobenius.norm(as.matrix(C.stat))^2
    
    if (method=="approximation") p.value <- 1-pchisq(Q.stat,p.ind1*p.ind2)
    if (method=="simulation" & scores=="rank") {Qs <-replicate(n.simu,.Q.simu.rank(rmvnorm(n,rep(0,p.ind1+p.ind2)),p1=p.ind1,p2=p.ind2,n=n));
                  p.value <- mean(Qs>Q.stat)}
    if (method=="simulation" & scores=="sign") {Qs <-replicate(n.simu,.Q.simu.sign(rmvnorm(n,rep(0,p.ind1+p.ind2)),p1=p.ind1,p2=p.ind2,n=n));
                  p.value <- mean(Qs>Q.stat)}
    if (method=="simulation" & scores=="normal") {Qs <-replicate(n.simu,.Q.simu.normal(rmvnorm(n,rep(0,p.ind1+p.ind2)),p1=p.ind1,p2=p.ind2,n=n));
                 p.value <-  mean(Qs>Q.stat)}
                 
    if (method=="permutation" & scores=="rank") {Qp <-replicate(n.simu,.Q.perm.rank(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n));
                  p.value <- mean(Qp>Q.stat)}
    if (method=="permutation" & scores=="sign") {Qp <-replicate(n.simu,.Q.perm.sign(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n));
                  p.value <- mean(Qp>Q.stat)}
    if (method=="permutation" & scores=="normal") {Qp <-replicate(n.simu,.Q.perm.normal(Z1.signs,Z1.ranks,Z2.signs,Z2.ranks,n));
                 p.value <-  mean(Qp>Q.stat)}
    STATISTIC<-Q.stat
    names(STATISTIC)<-"Q"
    
    Scores<-switch(scores, "sign"="Signs",
                  "rank"="Ranks",
                  "normal"="Normal Scores")  
    
    METHOD<-paste("Test of Independence Based on Marginal", Scores, "in an IC model")
    if(method=="approximation")
    {
    PARAMETER<-p.ind1*p.ind2
    names(PARAMETER)<-c("df")
    }
    else
    {
    PARAMETER<-n.simu
    names(PARAMETER)<-c("replications")
    }
    PVAL<-p.value
    res<-list(method=METHOD,statistic=STATISTIC,data.name=DNAME,parameter=PARAMETER,p.value=PVAL)
    class(res)<-"htest"
    return(res)
    
    }
