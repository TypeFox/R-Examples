summary.fregre.fd<-function(object,times.influ=3,times.sigma=3,draw=TRUE,...){
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
    isfdata<-is.fdata(y)
    n=nrow(x)
    if (!isfdata) {
     up=mean(object$residuals)+times.sigma*sqrt(object$sr2)
     lo=mean(object$residuals)-times.sigma*sqrt(object$sr2)
     i.atypical=which(object$residuals>up|object$residuals<lo)
     lim.influ=traza(object$H)/n
     influence=diag(object$H)
     i.influence=which(influence>times.influ*lim.influ)
     if (length(i.influence) == 0) i.influence=NA
     if (length(i.atypical) == 0) i.atypical=NA
     }
     if (object$call[[1]]=="fregre.pc") {
     if (object$lambda==0)     {
     cat(" *** Summary Functional Data Regression with Principal Components ***\n")
      object$lm$call<-object$call
      print(summary(object$lm))}
      else  {
     cat(" *** Summary Functional Ridge Regression with Principal Components*** \n\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
     cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
       cat("-Penalization parameter (lambda): ",object$lambda,"\n")
            }
            
#     object$lm$call<-object$call
#     print(summary(object$lm))
     var.1<-apply(object$fdata.comp$x, 2, var)
     pr.x= var.1/sum(var.1)
 cat("\n-With",length(object$l),"Principal Components is  explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n
-Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.ppc") {
     cat(" *** Summary Functional Regression with Penalized Principal Components ***\n")
      object$lm$call<-object$call
      print(summary(object$lm))
#            cat("\n-R squared: ",object$r2)
#      cat("\n-Residual variance: ",
#            object$sr2,"on ",n-object$df," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
       #     object$lm$call<-object$call
#     print(summary(object$lm))
     var.1<-apply(object$fdata.comp$x, 2, var)
     pr.x= var.1/sum(var.1)
 cat("\n-With",length(object$l),"Principal Components is explained ",round(sum(pr.x[object$l])*100
 ,2),"%\n of the variability of explicative variables. \n -Variability for each  principal components -PC- (%):\n")
    print(round(pr.x[object$l] * 100, 2))
    }
     if (object$call[[1]]=="fregre.pls") {
##     cat(" *** Summary Functional Data Regression with Partial Least Squares ***\n")
##      object$lm$call<-object$call
##      print(summary(object$lm))
     cat(" *** Summary Functional Regression with Partial Least Squares*** \n\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
              cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")

#     object$lm$call<-object$call
#     print(summary(object$lm))
#     var.1<-apply(object$fdata.comp$x, 2, var)
#     pr.x= var.1/sum(var.1)
# cat("\n-With",length(object$l),"Partial Least Squares is  explained ",round(sum(pr.x[object$l])*100
# ,2),"%\n of the variability of explicative variables. \n
#-Variability for each  Partial Least Squares -PLS- (%):\n")
#    print(round(pr.x[object$l] * 100, 2))
    }
    if (object$call[[1]]=="fregre.ppls") {
     cat(" *** Summary Functional Regression with Penalized Partial Least Squares ***\n")
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
              cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
       cat("-Lambda penalty: ",object$lambda)
      #     object$lm$call<-object$call
#     print(summary(object$lm))
#     var.1<-apply(object$fdata.comp$x, 2, var)
#     pr.x= var.1/sum(var.1)
# cat("\n-With",length(object$l),"Partial Least Squares is  explained ",round(sum(pr.x[object$l])*100
# ,2),"%\n of the variability of explicative variables. \n -Variability for each Partial Least Squares -PLS- (%):\n")
#    print(round(pr.x[object$l] * 100, 2))
    }

    if (object$call[[1]]=="fregre.basis") {
     cat(" *** Summary Functional Data Regression with representation in Basis *** \n")
     if (object$lambda==0)     {object$lm$call<-object$call
                                print(summary(object$lm))}
      else  {
            cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefs)
            cat("\n-R squared: ",object$r2)
            cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
            }
    }
  if (object$call[[1]]=="fregre.basis.cv") {
     cat(" *** Summary Functional Data Regression with representation in Basis *** \n\n")
      cat("-Call: ");    print(object$call)
            cat("\n")
            print(object$coefficients)
            cat("\n-R squared: ",object$r2)
#            cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
      cat("-Optimal Beta Basis: \n")
      print(object$basis.b.opt)
      cat("\n-Optimal lambda penalty=",object$lambda.opt,"\n")
#      print(object$Lfdobj)
    }
    if (object$call[[1]]=="fregre.np") {
     cat(" *** Summary Functional Non-linear Model *** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
#    cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
    }
  if (object$call[[1]]=="fregre.np.cv") {

     cat(" *** Summary Functional Non-linear Model *** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
#    cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
    }
     if (object$call[[1]]=="fregre.plm") {
     cat(" *** Summary Functional Semi-linear Model *** \n\n")
     cat("-Call: ");    print(object$call)
     cat("\n-Coefficients:  non functional covariates\n")
     print(object$coefficients)
      cat("\n-Bandwidth (h): ",object$h.opt)
    cat("\n-R squared: ",object$r2)
#    cat("\n-Residual variance: ",object$sr2,"\n")
cat("\n-Residual variance: ",
            object$sr2,"on ",n-object$df," degrees of freedom\n")
     }
    if (!isfdata) {
     cat("-Names of possible atypical curves: ");
     if (is.na(i.atypical[1]))     cat("No atypical curves \n")
     else   if (length(i.atypical)<11)  cat(rownames(x)[i.atypical],"\n")
           else cat(rownames(x)[i.atypical[1:10]],
           "\n It prints only the 10 most atypical curves. \n")
     cat("-Names of possible influence curves: ");
     if (is.na(i.influence[1])) cat("No influence curves \n")
     else  if (length(i.influence)<11) cat(rownames(x)[i.influence],"\n")
     else cat(rownames(x)[i.influence[1:10]],
     "\n It prints only the 10 most influence curves \n")
   }
   else  return(invisible(object)) #draw=FALSE
    if (draw) {
      oldpar <- par()
      C<-match.call()
      lenC=length(C)
      j=1
      while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1             }
        else {      j=j+1
                    ask=FALSE             }
       }
       if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
       else   par(mfrow=c(2,3))
 plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
     round(object$r2,2)))
 plot(object$fitted.values,object$residuals,ylab="Residuals",
    xlab="Fitted values",main="Residuals vs fitted.values")
    text(object$fitted.values[i.atypical],object$residuals[i.atypical],
    rownames(x)[i.atypical],cex=0.7)
    abline(h=mean(object$residuals),lwd=1,lty=2)
    abline(h=up,col=2,lwd=2,lty=2)
    abline(h=lo,col=2,lwd=2,lty=2)
#############
resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
main= "Scale-Location"
ylab23<-"Standardized residuals"
ylim <- c(0, max(resid.sd, na.rm = TRUE))
yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
plot(object$fitted.values,resid.sd, xlab = "Fitted values",
 ylab = yl, main = main,ylim = ylim)
 text(object$fitted.values[i.atypical],resid.sd[i.atypical],
 rownames(x)[i.atypical],cex=0.7)
 plot(diag(object$H),1:nrow(x),xlab="Leverage",ylab="Index.curves",
    main="Leverage")
text(diag(object$H)[i.influence],i.influence,
rownames(x)[i.influence],cex=0.7)
abline(v=times.influ*lim.influ,col=2,lwd=2,lty=2)
#  plot(density(object$residuals),main="Residuals")
    qqnorm(object$residuals,main="Residuals")
    boxplot(object$residuals,main="Residuals")
    par(mfrow=c(1,1))
    }
#    cat("\n")
return(invisible(list("Influence"=influence,"i.influence"=i.influence,
"i.atypical"=i.atypical)))
}
##############################################################################
##############################################################################
summary.fregre.gkam<-function(object,draw=TRUE,selec=NULL,times.influ=3,...){
   cat(" *** Summary Functional Data Regression with backfiting algorithm *** \n")
 print(object$family,"\n")
    nobs <- n<-length(object$y)
 cat("alpha=",as.numeric(format(object$alpha,digits=3)),"  n= ",n,"\n")
 cat("Algorithm converged?",ifelse(object$converged,"Yes","No")," Number  of iterations ",object$iter,"\n\n")
 res<- object$result
 namesres<-names(res)
 lenres<-length(res)
 eta<-object$linear.predictors
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
    if (object$family$family=="binomial") y<-as.numeric(factor(y,labels=c(0,1)))-1
    lenres<-length(object$result)
    tab<-matrix(NA,nrow=(lenres),ncol=3)
    rownames(tab)<-paste("f(",namesres,")",sep="")
    colnames(tab)<-c("h","cor(f(X),eta)","edf")
    cat("****    ****    ****    ****    ****    ****\n")
    for (i in 1:lenres) {
    hh<-res[[i]]$h.opt
    tab[i,1]<-as.numeric(format(hh,digits=3))
    tab[i,2]<-round(cor(object$effects[,i],object$fitted.values),3)
    tab[i,3]<-round(object$eqrank[[i]],1)
    }
    print(tab)
    cat("****    ****    ****    ****    ****    ****\n\n")
    cat("edf: Equivalent degrees of freedom\n")
    rowname<-rownames(res[[1]]$fdataobj)
    if (is.null(rowname)) rowname<-1:n
    else rowname<-rownames(x)
     H<-object$H
    influence=diag(H)
    lim.influ=sum(influence)/n
    i.influence=which(influence>times.influ*lim.influ)
    if (length(i.influence) == 0) i.influence=NA
    residual.df <- nobs - sum(object$eqrank)
    w<-object$prior.weights
    ycen = y - mean(y)
    r.sq <- round(1 - var(w * (y - object$fitted.values))*(nobs - 1)/(var(w * y)*residual.df),3)
    r.sq2 <- round( 1 - sum((object$fitted.values-y)^2)/sum(ycen^2),3)
    cat("Residual deviance=",round(object$deviance,3)," Null deviance=",round(object$null.deviance,3),"\n")
    dev.expl<-100*round((object$null.deviance - object$deviance)/object$null.deviance,3)
#     cat("R-sq. = ",1 - sum((object$residuals)^2)/sum(ycen^2),"  Deviance explained = ", dev.expl,"%","\n")
    cat("AIC= ",round(object$aic,3),"Deviance explained=",dev.expl,"%","\n")
    cat("R-sq.=",r.sq2," R-sq.(adj)=",r.sq,"\n")
    cat("Names of possible influence curves: ");
    if (is.na(i.influence[1])) cat("No influence curves \n")
    else  if (length(i.influence)<11) cat(rowname[i.influence],"\n")
    else cat(rowname[i.influence[1:10]],
     "\n It prints only the 10 most influence curves \n")
    cat("\n")
    if (draw) {
      le<-lenres
      C<-match.call()
      lenC=length(C)
      j=1
      while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1             }
        else {      j=j+1
                    ask=FALSE             }
       }

          ############
       if (ask) {          par(mfrow=c(1,1))
          }
       else                par(mfrow = c(2,3))


# plot(object$fitted.values,y,xlab="Fitted values",main=paste("R-squared=",
#     round(r.sq,2)))
col1<-rep(1,n)
#if (is.factor(y)) col1=as.numeric(y)
#    if (family$family=="binomial" & !is.factor(y)) y<-factor(y)
    if (object$family$family=="binomial") col1=y+2

    plot(object$linear.predictors,y,xlab="Linear predictors",main=paste("R-sq=",
     round(r.sq2,2)))
    points(object$linear.predictors,object$fitted.values,col=col1)
    legend(min(object$linear.predictors),0.95,legend=c("y","fitted values"),col=1:2,lty=1:2,
    box.col="white",xjust=0)
    plot(object$fitted.values,object$residuals,ylab="Residuals",
    xlab="Fitted values",main="Residuals vs fitted.values",col=col1)
    abline(h=mean(object$residuals),lwd=1,lty=2)
#############
resid.sd=sqrt(abs(object$residuals/sd(object$residuals)))
main= "Scale-Location"
ylab23<-"Standardized residuals"
ylim <- c(0, max(resid.sd, na.rm = TRUE))
yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
plot(object$fitted.values,resid.sd, xlab = "Fitted values",
ylab = yl, main = main,ylim = ylim,col=col1)
# text(object$fitted.values[i.atypical],resid.sd[i.atypical],rowname[i.atypical],cex=0.7)
 plot(diag(H),1:n,xlab="Leverage",ylab="Index.curves",main="Leverage",col=col1)
text(influence[i.influence],i.influence,rowname[i.influence],cex=0.7)
 abline(v=times.influ*lim.influ,col=4,lwd=2,lty=2)
    qqnorm(object$residuals,main="Residuals",col=col1)
    boxplot(object$residuals,main="Residuals")
    X<-object$effects
    if (length(table(y) > 6)) {colores = 1}
    else { colores = y + 1 }
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
    if (!ask) {
              dev.new()
              cc<-if(le>6) {cc<-c(ceiling(le/3), 3)}
              else {
              if (le>=2) cc<-c(ceiling(le/2), 2)
              else cc<-c(1,1)
              }
              par(mfrow = cc)
    }
    else {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
    if (!is.null(selec)) {
                par(mfrow=c(1,1))
                plot(X[, selec], eta, col = col1, ylab = "Linear predictors",
              xlab = paste("f(", namesres[selec], ")", sep = ""),
              main = paste(namesres[selec], "edf:", round(object$eqrank[selec],1)))
               if (length(table(y)) == 2) {
                abline(h = 0); abline(v = 0)
                    }
    }
    else {
      for (i in 1:(ncol(X)-1)) {
         plot(X[, i], eta, col = col1, ylab = "Linear predictors",
              xlab = paste("f(", namesres[i], ")", sep = ""),
                main = paste(namesres[i], "edf:", round(object$eqrank[i],1)))
               if (length(table(y)) == 2) {
                abline(h = 0);abline(v = 0)
                }
                }
    }
    par(mfrow=c(1,1))
    }
    cat("\n")
return(invisible(list("Influence"=influence,"object"=object)))
}

##############################################################################
##############################################################################
print.fregre.gkam<-function(x,digits = max(3, getOption("digits") - 3),...){
 object<-x
 print(object$family,"\n")
# cat("Algorithm converged?",ifelse(object$converged,"Yes","No")," Number  of iterations ",object$iter,"\n")
 res<- object$result
 namesres<-names(res)
 lenres<-length(res)
    x<-object$fdataobj$data
    t=object$fdataobj$argvals
    y<-object$y
   if (object$family$family=="binomial") y<-as.numeric(factor(y,labels=c(0,1)))-1
    
    nobs <- n<-length(y)
    lenres<-length(object$result)
    tab<-matrix(NA,nrow=lenres,ncol=3)
    rownames(tab)<-paste("f(",namesres,")",sep="")
    colnames(tab)<-c("h","cor(f(X),eta)","edf")
    cat("alpha=",as.numeric(format(object$alpha,digits=3)),"  n= ",n,"\n")
    cat("****    ****    ****    ****    ****    ****\n")
    for (i in 1:lenres) {
    hh<-res[[i]]$h.opt
    tab[i,1]<-as.numeric(format(hh,digits=3))
    tab[i,2]<-round(cor(object$effects[,i],object$fitted.values),3)
    tab[i,3]<-round(object$eqrank[i],1)
    }
    print(tab)
    cat("****    ****    ****    ****    ****    ****\n")
    cat("edf: Equivalent degrees of freedom\n")
    residual.df <- nobs - sum(object$eqrank)
    w<-object$prior.weights
    ycen = y - mean(y)
    r.sq <- round(1 - var(w * (y - object$fitted.values))*(nobs - 1)/(var(w * y)*residual.df),3)
    r.sq2 <- round( 1 - sum((object$fitted.values-y)^2)/sum(ycen^2),3)
#    cat("Residual deviance= ",round(object$deviance,3),"  Null deviance= ",round(object$null.deviance,3),"\n")
    dev.expl<-100*round((object$null.deviance - object$deviance)/object$null.deviance,3)
#     cat("R-sq. = ",1 - sum((object$residuals)^2)/sum(ycen^2),"  Deviance explained = ", dev.expl,"%","\n")
    cat("AIC=",round(object$aic,3),"  Deviance explained =", dev.expl,"%","\n")
    cat("R-sq.=",r.sq2,"  R-sq.(adj)=",r.sq)
    cat("\n")
return(invisible(object))
}

##############################################################################
##############################################################################
kgam.H<-function(object,inverse="svd") {
#print("kgam.H")
 lenH<-length(object)
 if (lenH==1) return(object[[1]]$H)
 else {
 SS.list<-SS2<-list()
 n<-ncol(object[[1]]$H)
 SS<-matrix(NA,ncol=(lenH)*n,nrow=(lenH)*n)
 II=diag(n)
 unos<-matrix(1,ncol=n,nrow=n)/n
 M<-object[[1]]$H
# MM=sweep(M,1,apply(M,1,mean),"-")
 MM=sweep(M,1,rowMeans(M),"-")
 if (lenH>1) {
  for (i in 2:lenH) {
#   MMaux=sweep(object[[i]]$H,1,apply(object[[i]]$H,1,mean),"-")
  MMaux=sweep(object[[i]]$H,1,rowMeans(object[[i]]$H),"-")
   MM<-rbind(MM,MMaux)
   }
  }
  MM1<-matrix(rep(MM,lenH),nrow=(n*lenH))
  DD<-kronecker(diag(lenH),outer(rep(1,n),rep(1,n)))
  D1<-abs(DD-1)
  SS<-MM1*D1+diag(n*lenH)
#  cat("kappa=",kappa(SS),"\n")    
    slv <- try(solve(SS),silent=TRUE)
#    print(slv)
    if (class(slv)=="try-error") {
     sv<-svd(SS)
     slv<-drop((sv$v%*%diag(1/sv$d)%*%t(sv$u)))
    warning("Inverse of sigma computed by SVD")
    }    
  SSinv<-switch (inverse,solve=slv,
               svd={
                       res.X=svd(SS)
                      (res.X$v)%*%diag(res.X$d^(-1))%*%t(res.X$u)})
 H2=SSinv%*%MM
 HH<-unos+H2[1:n,1:n]
 if (lenH>1) {
  for (i in 2:lenH) {HH=HH+H2[(n*(i-1)+1):(n*i),]  }
  }
  HH
  }
 }
##############################################################################
##############################################################################
 