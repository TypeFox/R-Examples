summary.lmmlasso <- function(object,...)
 {
    if (object$stopped)
      {
        cat("Since the algorithm stopped due to","\n")
        cat("|active.set|>=min(p,ntot)","\n")
        cat("no summary is available.","\n")
        cat("Set stopSaturated=FALSE or increase lambda!","\n")
      } else {
    # 1. part
    table1 <- round(c(object$aic,object$bic,object$logLik,object$deviance,object$objective),1)
    names(table1) <- c("AIC","BIC","logLik","deviance","objective")
    cat("Model fitted by",object$method,"for","lambda =",round(object$lambda,3),":\n")
    print(table1)

    # 2a. part
    cat("","\n")
    cat("Random effects:",object$pdMat,"\n")
    tableVar <- c(diag(object$Psi),object$sigma^2)
    tableSd <- c(sqrt(diag(object$Psi)),object$sigma)
    matrixVar <- cbind(tableVar,tableSd)
    colnames(matrixVar) <- c("Variance","Std.Dev.")

    znames <- colnames(object$data$z)

    if (is.null(znames))
      {
        if (length(object$ranInd)>1) rownames(matrixVar) <- c(paste("X",object$ranInd,sep=""),"Residual")
        else rownames(matrixVar) <-c(paste("X",object$ranInd,sep=""),"Residual")
      } else
      {
        rownames(matrixVar) <- c(colnames(object$data$z),"Residual")
      }
    
    print(matrixVar)

    if ((object$pdMat=="pdSym")&(length(object$pars)>1))
      {
        cat("","\n")
        cat("Random effects: Correlations","\n")
        matrixCor <- object$corPsi
        if (dim(object$data$z)[[2]]>1) colnames(matrixCor) <- rownames(matrixCor) <- c("(Intercept)",paste("X", 2:dim(object$data$z)[[2]],sep=""))
        print(matrixCor)
      }
    
    # 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    penalty <- rep("",length(object$coefficients))
    penalty[object$nonpen] <- "(n)"
    table3 <- data.frame(object$coefficients,penalty)

    xnames <- colnames(object$data$x)

    if (is.null(xnames))
     {
       rownames(table3) <- c("(Intercept)",paste("X", 2:length(object$coefficients),sep=""))
     } else
     {
       rownames(table3) <- c("(Intercept)",colnames(object$data$x)[-1]) 	
     }
    table3 <- table3[object$coefficients!=0,]
    colnames(table3) <- c("Estimate"," ")
    cat("|active set|=",sum(object$coefficients!=0),"\n")
    print(table3)

    # 4th part
    cat("","\n")
    cat("Number of iterations:",object$counter,"\n")
  }
 }

print.lmmlasso <- function(x,...) {
    # 0. part
    print(x$call)
  
    # 2a. part
    cat("","\n")
    cat("Random effects:",x$pdMat,"\n")
    tableVar <- c(diag(x$Psi),x$sigma^2)
    tableSd <- c(sqrt(diag(x$Psi)),x$sigma)
    matrixVar <- cbind(tableVar,tableSd)
    colnames(matrixVar) <- c("Variance","Std.Dev.")
    if (dim(x$data$z)[[2]]>1) rownames(matrixVar) <- c("(Intercept)",paste("X", 2:dim(x$data$z)[[2]],sep=""),"Residual")
    else rownames(matrixVar) <-c("(Intercept)","Residual")
    print(matrixVar)
    
    # 3. part
    cat("","\n")
    cat("Fixed effects:","\n")
    cat("|active set|=",sum(x$coefficients!=0),"\n")

}

plot.lmmlasso <- function(x,...) {

par(mfrow=c(3,2))

# Tukey-Anscombe plot
plot(x$residuals~x$fitted.values,col=x$data$grp,main="Tukey-Anscombe Plot",xlab="fitted values",ylab="raw residuals")
abline(h=0,col="grey")

# QQ-Plot of the residuals
qqnorm(x$residuals,col=x$data$grp,main="QQ-Plot of the residuals")
qqline(x$residuals)

# QQ-plot of the random effects
ranef <- unlist(x$random,recursive=FALSE)
ranef.sd <- rep(sqrt(diag(x$Psi)),length(levels(x$data$grp)))
ranef.st <- ranef/ranef.sd
qqnorm(ranef.st,main="QQ-Plot of the standardized random effects",col=rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp))))
qqline(ranef.st)

# boxplots of the random effects
boxplot(ranef~as.factor(rep(1:dim(x$data$z)[[2]],length(levels(x$data$grp)))),border=1:dim(x$data$z)[[2]],
        main="Random effects by effects",xlab="random effects index")
abline(h=0,col="grey",lty=2)

# fitted vs. observed y
plot(x$fitted,x$data$y,xlab="fitted values",ylab="observed values",main="response ~ fitted")
abline(a=0,b=1)

# histogram of the fixed effects
hist(x$coefficients,main=paste("|active set| = ",sum(x$coefficients!=0)),xlab="fixed effects")
rug(x$coefficients,lwd=3)

}
