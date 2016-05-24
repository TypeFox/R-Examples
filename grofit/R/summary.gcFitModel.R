summary.gcFitModel <-
function(object, ...)
{
# object of class gcFitModel

contents.fitted.param     = c("mu.model", "lambda.model", "A.model", "integral.model",
                               "stdmu.model", "stdlambda.model", "stdA.model",
                               "ci90.mu.model.lo", "ci90.mu.model.up",
                               "ci90.lambda.model.lo", "ci90.lambda.model.up",
                               "ci90.A.model.lo", "ci90.A.model.up",
                               "ci95.mu.model.lo", "ci95.mu.model.up",
                               "ci95.lambda.model.lo", "ci95.lambda.model.up",
                               "ci95.A.model.lo", "ci95.A.model.up")


if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){


   table<-rep(NA,length(contents.fitted.param))
}
else{
   table <- c(object$parameter$mu[1], object$parameter$lambda[1],  object$parameter$A[1], object$parameter$integral,
           object$parameter$mu[2], object$parameter$lambda[2],  object$parameter$A[2], 
           object$parameter$mu[1]-1.645*object$parameter$mu[2], object$parameter$mu[1]+1.645*object$parameter$mu[2],
           object$parameter$lambda[1]-1.645*object$parameter$lambda[2], object$parameter$lambda[1]+1.645*object$parameter$lambda[2],
           object$parameter$A[1]-1.645*object$parameter$A[2], object$parameter$A[1]+1.645*object$parameter$A[2],
           object$parameter$mu[1]-1.96*object$parameter$mu[2], object$parameter$mu[1]+1.96*object$parameter$mu[2],
           object$parameter$lambda[1]-1.96*object$parameter$lambda[2], object$parameter$lambda[1]+1.96*object$parameter$lambda[2],
           object$parameter$A[1]-1.96*object$parameter$A[2], object$parameter$A[1]+1.96*object$parameter$A[2])

}
table <- data.frame(t(table))
colnames(table) <- contents.fitted.param
summary.gcFitModel <- table

}

