summary.gcBootSpline <-
function(object, ...)
{
# object of class gcBootSpline
contents.bootstrap        <- c("mu.bt", "lambda.bt", "A.bt", "integral.bt", "stdmu.bt", "stdlambda.bt", "stdA.bt", "stdintegral.bt",
                               "ci90.mu.bt.lo", "ci90.mu.bt.up", "ci90.lambda.bt.lo", "ci90.lambda.bt.up",
                               "ci90.A.bt.lo", "ci90.A.bt.up", "ci90.integral.bt.lo", "ci90.integral.bt.up",
                               "ci95.mu.bt.lo", "ci95.mu.bt.up", "ci95.lambda.bt.lo", "ci95.lambda.bt.up",
                               "ci95.A.bt.lo", "ci95.A.bt.up", "ci95.integral.bt.lo", "ci95.integral.bt.up")

							   
if (object$bootFlag==FALSE){
    table<-rep(NA,length(contents.bootstrap))
}
else{
     mu          <- mean(object$mu, na.rm=TRUE)
     lambda      <- mean(object$lambda, na.rm=TRUE)
     A           <- mean(object$A, na.rm=TRUE)
     integral    <- mean(object$integral, na.rm=TRUE)

     mu.sd       <- sd(object$mu, na.rm=TRUE)
     lambda.sd   <- sd(object$lambda, na.rm=TRUE)
     A.sd        <- sd(object$A, na.rm=TRUE)
     integral.sd <- sd(object$integral, na.rm=TRUE)

     table <- c(mu, lambda, A, integral, mu.sd, lambda.sd, A.sd, integral.sd,
                mu-1.645*mu.sd, mu+1.645*mu.sd, lambda-1.645*lambda.sd,     lambda+1.645*lambda.sd,
                A-1.645*A.sd,   A+1.645*A.sd,   integral-1.645*integral.sd, integral+1.645*integral.sd,
                mu-1.96*mu.sd,  mu+1.96*mu.sd,  lambda-1.96*lambda.sd,      lambda+1.96*lambda.sd,
                A-1.96*A.sd,    A+1.96*A.sd,    integral-1.96*integral.sd,  integral+1.96*integral.sd)
}

table               <- data.frame(t(table))
colnames(table)     <- contents.bootstrap
summary.gcBootSpline <- table

}

