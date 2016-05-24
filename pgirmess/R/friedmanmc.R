friedmanmc<-function(y,groups,blocks,probs=0.05){

    if ((is.matrix(y)) | (is.data.frame(y))){
        groups <- factor(c(col(y)))
        blocks <- factor(c(row(y)))
    }                          
    
    else {
        if (any(is.na(groups)) || any(is.na(blocks))) 
            stop("NA's are not allowed in groups or blocks")
        if (any(diff(c(length(y), length(groups), length(blocks))))) 
            stop("y, groups and blocks must have the same length")
        if (any(table(groups, blocks) != 1)) 
            stop("Not an unreplicated complete block design")
        groups <- factor(groups)
        blocks <- factor(blocks)
    }

    N<-length(levels(blocks));k<-length(levels(groups))
    y <- y[order(groups,blocks)]
    mat<-matrix(y,nrow=N,ncol=k,byrow=FALSE)
    for (i in 1:length(mat[,1])) mat[i,]<-rank(mat[i,])
    sumrank<-colSums(mat)
    names(sumrank)<-levels(groups)
    dif<-NULL; noms<-NULL
    for (i in 1:(length(sumrank)-1)){
        for (j in (i+1):length(sumrank)){            
            dif<-c(dif,abs(sumrank[i]-sumrank[j]))
            noms<-c(noms,paste(names(sumrank[i]),"-",names(sumrank[j]),sep=""))
        }
    }

    names(dif)<-noms
    dif1<-dif
    alpha<-qnorm(probs/(k*(k-1)),lower.tail=FALSE)*sqrt((N*k*(k+1))/6)
    meth<-"Multiple comparisons between groups after Friedman test"
    output<-list(statistic=meth, p.value=probs, dif.com=data.frame(obs.dif=dif1,critical.dif=alpha,difference=ifelse((dif1-alpha)>0,TRUE,FALSE)))
    class(output)<-c("mc","list")
    output
}
