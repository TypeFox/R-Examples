`relabel` <-
function(cls, print.loss=TRUE){
    
     min.perm.labels <- function(cl,Q){ 
        cost.mat <- apply(matrix(1:nclus),1,function(j) apply(Q,2, function(x) sum(x[cl==j])))
        sol.assign <- lp.transport(cost.mat,row.signs=rep("=",nclus), row.rhs=rep(1,nclus), col.signs=rep("=",nclus),col.rhs=rep(1,nclus))
        perm <- unlist(apply(sol.assign$solution,2,function(x) which(x > 0.1)))
        list(cl=perm[cl], val= sol.assign$objval)
    }
    
    k <- apply(cls,1,function(x) length(table(x)))
    if(length(table(k))>1){
        k.max <- as.numeric(names(table(k))[which.max(table(k))])
        cls <- cls[k==k.max, ]
        warning("only clusterings with the most common number of groups used for relabelling")
    }
    
   
    nclus <- length(table(cls))
    n <- ncol(cls)
    M <- nrow(cls)
    
    loss.val <- log(nclus)*n*M  # Gibt Verlustfunktionswert, wenn Gleichverteilung für jede Beobachtung vorliegt
    if(print.loss) cat("Value Loss Function:", loss.val,"\n")
    repeat{
        # Step 1
        P <- t(apply(cls, 2, function(x) tabulate(x,nbins=nclus)))/M 
        Q <- -log(P)
        Q[Q==Inf] <- log(100*M) # lp.transport kann nicht mit Inf umgehen
        # Step2
        res.perm <- t(apply(cls, 1, min.perm.labels,Q=Q))
        cls.new <- matrix(unlist(lapply(res.perm, function(x)x$cl)), ncol=n,byrow=TRUE)
        vals <- unlist(lapply(res.perm, function(x) x$val))
        loss.val.new <- sum(vals)
        cls <- cls.new
        if(print.loss) cat("Value Loss Function:", loss.val.new,"\n")
        if(abs(loss.val.new -loss.val)< 1E-6) break
        loss.val <- loss.val.new
    }
list(cls=cls, P=P, loss.val=loss.val, cl=apply(P,1,which.max))
}
