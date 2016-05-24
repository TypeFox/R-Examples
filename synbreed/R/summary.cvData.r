# summary method for cross validation results
summary.cvData <- function(object,...){

     obj <- object
     ans <- list()

     # summary of method
     ans$k <- obj$k
     ans$Rep <- obj$Rep
     ans$est.method <- obj$VC.est.method
     ans$sampling <- obj$sampling
     if(obj$sampling=="commit") ans$sampling <- "committed"
     ans$nr.ranEff <- obj$nr.ranEff

     # number of individuals
     if(!is.null(obj$n.DS)){
     ans$nmin.DS <- min(obj$n.DS)
     ans$nmax.DS <- max(obj$n.DS)
     }

     # size of TS
     ans$nmin.TS <- min(obj$n.TS)
     ans$nmax.TS <- max(obj$n.TS)

     # Results
     # Predictive ability
     colmean.pa <- colMeans(obj$PredAbi)
     ans$se.pa <- format(sd(colmean.pa)/sqrt(length(colmean.pa)),digits=4,nsmall=4)
     ans$min.pa <- format(min(obj$PredAbi),digits=4,nsmall=4)
     ans$mean.pa <- format(mean(obj$PredAbi),digits=4,nsmall=4)
     ans$max.pa <- format(max(obj$PredAbi),digits=4,nsmall=4)
     # Spearman's rank correlation 
     if(!is.null(obj$rankCor)){
     colmean.rc <- colMeans(obj$rankCor)
     ans$se.rc <- format(sd(colmean.rc)/sqrt(length(colmean.rc)),digits=4,nsmall=4)
     ans$min.rc <- format(min(obj$rankCor),digits=4,nsmall=4)
     ans$mean.rc <- format(mean(obj$rankCor),digits=4,nsmall=4)
     ans$max.rc <- format(max(obj$rankCor),digits=4,nsmall=4)
     }
     # Bias
     colmean.b <- colMeans(obj$bias)
     ans$se.b <- format(sd(colmean.b)/sqrt(length(colmean.b)),digits=4,nsmall=4)
     ans$min.b<- format(min(obj$bias),digits=4,nsmall=4)
     ans$mean.b <- format(mean(obj$bias),digits=4,nsmall=4)
     ans$max.b <-format(max(obj$bias),nsmall=4,digits=4)
     # 10% best
     if(!is.null(obj$m10)){
     ans$se.10 <- format(sd(as.vector(obj$m10))/sqrt(length(obj$m10)),digits=4,nsmall=4)
     ans$min.10<- format(min(obj$m10),digits=2,nsmall=2)
     ans$mean.10 <- format(mean(obj$m10),digits=2,nsmall=2)
     ans$max.10 <-format(max(obj$m10),nsmall=2,digits=2)
     }
     # Mean squared error
     if(!is.null(obj$mse)){
     colmean.mse <- colMeans(obj$mse)
     ans$se.mse <- format(sd(colmean.mse)/sqrt(length(colmean.mse)),digits=4,nsmall=4)
     ans$min.mse <- format(min(obj$mse),digits=3,nsmall=3)
     ans$mean.mse <- format(mean(obj$mse),digits=3,nsmall=3)
     ans$max.mse <- format(max(obj$mse),digits=3,nsmall=3)
     }
     # Seed
     ans$Seed <- obj$Seed
     ans$rep.seed <- obj$rep.seed

     class(ans) <- "summary.cvData"
     ans
}

# print method for summary.pedigree
print.summary.cvData <- function(x,...){
    cat("Object of class 'cvData' \n")
    cat("\n",x$k,"-fold cross validation with",x$Rep,"replication(s) \n")
    cat("     Sampling:                ",x$sampling,"\n")
    cat("     Variance components:     ",x$est.method,"\n")
    cat("     Number of random effects:",x$nr.ranEff,"\n")
    if(!is.null(x$nmin.DS)) cat("     Number of individuals:   ",x$nmin.DS,"--", x$nmax.DS," \n")
    cat("     Size of the TS:          ",x$nmin.TS,"--", x$nmax.TS," \n")
    cat("\nResults: \n")
    cat("                      Min \t  Mean +- pooled SE \t  Max \n")
    cat(" Predictive ability: ",x$min.pa," \t ",x$mean.pa,"+-",x$se.pa," \t ",x$max.pa,"\n")
    if(!is.null(x$mean.rc)) cat(" Rank correlation:   ",x$min.rc," \t ",x$mean.rc,"+-",x$se.rc," \t ",x$max.rc,"\n")
    if(!is.null(x$mean.mse)) cat(" Mean squared error: ",x$min.mse," \t ",x$mean.mse,"+-",x$se.mse," \t ",x$max.mse,"\n")
    cat(" Bias (reg. slope)   ",x$min.b," \t ",x$mean.b,"+-",x$se.b," \t ",x$max.b,"\n")
    if(!is.null(x$mean.10)) cat(" 10% best predicted: ",x$min.10," \t ",x$mean.10,"+-",x$se.10," \t ",x$max.10,"\n")
    cat("\nSeed start: ",x$Seed,"\n")
    cat("Seed replications: \n")
    print(x$rep.seed)
} 



