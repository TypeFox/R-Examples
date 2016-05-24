print.GPvam <-
function (x, ...) 
{
    cat("Object of class 'GPvam'. Use functions 'plot' and 'summary'.", "\n") 
    cat("Object contains elements:  loglik, teach.effects, parameters Hessian, R_i, teach.cov, mresid, cresid, y, yhat, stu.cov, num.obs, num.student, num.year, num.teach, yhat.m, sresid, yhat.s\n")
    cat("Number of iterations: ",x$iter,"\n",sep="")
    cat("Log-likelihood: ",x$loglik,"\n",sep="")
    cat("Parameter estimates:\n")
    print(x$parameters)
    cat("Teacher Effects\n")
    print(x$teach.effects)
    cat("\n")
    
}
