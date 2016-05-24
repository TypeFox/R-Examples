summary.semsfa<-function(object,...){
    cat("\nSemiparametric Stochastic Frontier Model\n")    
    print(object$call)
    e.coef<-c(object$lambda,object$sigma)
    if(object$n.boot>0) 
    {
     b.se<-object$b.se    
     b.t<-e.coef/b.se
     b.pv<-2 * pt(abs(b.t), df = object$residual.df, lower.tail = FALSE)
    }else 
    {
     b.se<-b.t<-b.pv<-rep(NA,2)    
    }
    cat("\nCONDITIONAL EXPECTATION ESTIMATE\n")
    print(summary(object$reg))
    cat("\nVARIANCE COMPONENTS ESTIMATE (boostrap replicates =",object$n.boot,")")
    p.table <- cbind(e.coef, b.se, b.t, b.pv)
    dimnames(p.table) <- list(c("lambda", "sigma"), c("Estimate", 
        "Boot SE", "t value", "Pr(>|t|)"))
    object <- list(object,coefficients=p.table)#call=object$call,
    cat("\n")
    printCoefmat(object$coefficients, P.values=TRUE,has.Pvalue=TRUE)
    object$b.t<- b.t
    object$b.pv<-b.pv
    class(object)<-"summary.semsfa"
    return(object)    
}   