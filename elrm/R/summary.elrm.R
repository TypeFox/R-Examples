`summary.elrm` <-
function(object,...)
{
    inferences = as.data.frame(cbind(round(as.numeric(object$coeffs),5),round(as.numeric(object$p.values),5),round(as.numeric(object$p.values.se),5),object$mc.size));
    
    results = data.frame(row.names=names(object$coeffs), inferences);
    names(results) = c("estimate","p-value","p-value_se","mc_size");
    
    cat("\nCall:\n");
    print(object$call.history);
    cat('\n');
    cat("Results:\n");
    print(results,quote=F);
    cat('\n');
    cat(object$ci.level,"% Confidence Intervals for Parameters\n",sep="");
    cat('\n');
    print(object$coeffs.ci);
    cat('\n');
}

