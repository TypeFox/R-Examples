#' Generalised linear models with variable selection
#'
#' \code{glmsel} is used to fit generalised linear models with optionally performing variable selection using stepwise regression, lasso or elastic net methods.
#' @usage glmsel(formula, data=environment(), varsel=FALSE, criterion="AIC",
#' direction="backward", indices=NULL, train=0.3, family, enet.alpha=0.5)
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment containg the variables in the model. If not specified, the variables are taken from the current environment.
#' @param varsel a method of variable selection to be used. The default is \code{"FALSE"}. Available methods include: stepwise regression \code{"step"}, LASSO \code{"lasso"}, elastic net\code{"enet"}.
#' @param criterion when \code{varsel="step", criterion} allows to select a method of calculating statistic for model comparison. The default is  \code{"AIC"}. Less liberal, BIC penalty can be used by typing \code{"BIC"}.
#' @param direction the mode of stepwise search, can be one of \code{"both"}, \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}. If the scope argument is missing the default for \code{direction} is \code{"backward"}.
#' @param indices vector of \code{0} and \code{1} values indicating which observations are to be used as train and test when \code{varsel="lasso" or "enet"}.
#' @param train if \code{indices=NULL}, the function will randomly assign observations as train and test. \code{train} specifies what percentage of data will be used as train observations. Can take values from \code{0.1} to \code{0.9}.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function.
#' @param enet.alpha The elastic net mixing parameter, with 0=a= 1. The penalty is defined as\deqn{(1-a)/2||\beta||\_{2}^2+a||\beta||_1}\code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty. The default value is \code{0.5}.
#' @return A "glmsel" object is returned, for which print, plot and summary methods can be used.
#' @examples data(bank)
#' glmsel(y~age+as.factor(job)+as.factor(marital)+as.factor(education)+as.factor(default)+balance
#' +as.factor(housing)+as.factor(loan)+as.factor(contact)+as.factor(day)+month+duration+campaign
#' +pdays+previous, data=bank, varsel="enet")
#' @import stats
#' @import glmnet
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
glmsel<-function(formula, data=environment(), varsel=FALSE, criterion="AIC", direction="backward", indices=NULL, train=0.3, family="gaussian",
                 enet.alpha=0.5) UseMethod("glmsel")

#' @export
glmsel.default <-function(formula, data=environment(), varsel=FALSE, criterion="AIC", direction="backward", indices=NULL, train=0.3, family="gaussian",
                          enet.alpha=0.5)
{
  if (!missing(data)) data<- as.data.frame(data)
  x<- model.matrix(formula, data=data)
  y<-model.frame(formula, data=data)[,1]
  if(!is.numeric(y)) stop("response vector is not numeric")
  if(varsel=="enet" & enet.alpha==0) warning("enet.alpha=0 will perform the lasso fit")
  if(train>=1 | train<=0) stop("train must be > 0 and < 1")
  if(enet.alpha<0) stop("enet.alpha cannot be negative")
  if(!all(indices %in% c(0,1))) stop("indices column should only have values 0 and 1")
  if(is.null(indices)==FALSE & all(indices %in% 0)) stop("indices column cannot contain only 0 values")
  if(is.null(indices)==FALSE & all(indices %in% 1)) stop("indices column cannot contain only 1 values")
  if(is.null(indices) & (!missing(varsel))){indices<-numeric(nrow(x)) #obtain indices for data splitting
  if((family=="exponential" | family=="Exponential") & (varsel=="enet" | varsel=="lasso")) stop("lasso and enet variable selection methods don't work with exponential family")
  r<-sample(1:nrow(x), size=train*nrow(x))
  indices[r]<-1
  indices[-r]<-0
  }#obtain indices for data splitting
  var.rejected <- character(0)
  #VARIABLE SELECTION
  if(varsel=="lasso"){
    if(colnames(x)[1]=="(Intercept)"){#here we do need to exclude intercept
      x.test<- x[indices==0,-1]
      x.train<- x[indices==1,-1]
      y.train<- y[indices==1]
      intercept<-TRUE
    }else{x.test<- x[indices==0,]
    x.train<- x[indices==1,]
    y.train<- y[indices==1]
    intercept<-FALSE
    }
    cv <- cv.glmnet(x=x.train, y=y.train, family=family, alpha=1)
    gn<-glmnet(x=x.train, y=y.train, family=family, alpha=1)
    var.all <- predict.cv.glmnet(cv, type="coefficients")[,1]#here x should be test.x #rerun this line with xnew=x.test and type not equal to coefficients to get the estimates for regression
    var.all<-var.all[names(var.all)!="(Intercept)"]
    #handle selected factors
    var<-names(var.all)[which(var.all!=0)]
    n<- grep( "as[.]factor.*", var)
    var[n] <- gsub( ").+$", ")", var[n] )
    var<- unique(var)
    #variables rejected
    var.all2<-names(var.all)
    n2<- grep( "as[.]factor.*", var.all2)
    var.all2[n2] <- gsub( ").+$", ")", var.all2[n2])
    var.all2<- unique(var.all2)
    var.rejected<-setdiff(var.all2, var)
    #transform vectors with variables to a formula object
    if(length(var)==0 & intercept==FALSE) stop("All variables were rejected by the algorithm")#case when all variables are rejected and no intercept was chosen
    if(length(var)==0) var<-"1"#case when all variables are rejected
    if(intercept==TRUE){formula<-formula(paste(formula[[2]], "~",paste(var, collapse=" + ")))
    }else{formula<-formula(paste(formula[[2]], "~",paste(var, collapse=" + "),"-1"))
    }
  }
  if(varsel=="enet"){
    if(colnames(x)[1]=="(Intercept)"){#create case when intercept is not present
      x.test<- x[indices==0,-1]
      x.train<- x[indices==1,-1]
      y.train<- y[indices==1]
      intercept<-TRUE
    }else{x.test<- x[indices==0,]
    x.train<- x[indices==1,]
    y.train<- y[indices==1]
    intercept<-FALSE
    }
    cv <- cv.glmnet(x=x.train, y=y.train, family=family, alpha=enet.alpha)
    gn<-glmnet(x=x.train, y=y.train, family=family, alpha=enet.alpha)
    var.all <- predict.cv.glmnet(cv, type="coefficients")[,1]#here x should be test.x #rerun this line with xnew=x.test and type not equal to coefficients to get the estimates for regression
    var.all<-var.all[names(var.all)!="(Intercept)"]
    #handle selected factors
    var<-names(var.all)[which(var.all!=0)]
    n<- grep( "as[.]factor.*", var)
    var[n] <- gsub( ").+$", ")", var[n] )
    var<- unique(var)
    #variables rejected
    var.all2<-names(var.all)
    n2<- grep( "as[.]factor.*", var.all2)
    var.all2[n2] <- gsub( ").+$", ")", var.all2[n2])
    var.all2<- unique(var.all2)
    var.rejected<-setdiff(var.all2, var)
    #transform vectors with variables to a formula object
    if(length(var)==0 & intercept==FALSE) stop("All variables were rejected by the algorithm")#case when all variables are rejected and no intercept was chosen
    if(length(var)==0) var<-"1"#case when all variables are rejected
    if(intercept==TRUE){formula<-formula(paste(formula[[2]], "~",paste(var, collapse=" + ")))
    }else{formula<-formula(paste(formula[[2]], "~",paste(var, collapse=" + "),"-1"))
    }
  }
  if(varsel=="step"){
    if(criterion=="BIC") k<-log(nrow(x))
    if(criterion=="AIC") k<-2
    myglm<-step(glm(formula=formula, family=family, data=data), k=k, trace=FALSE)
    var.rejected<-gsub("^- ", "",(myglm$anova$Step)[-1])
    #get var rejected thingy
    #if(length(var.rejected)==0){myglm$var.rejected<-"no variables were rejected"
  }else{
    myglm<-glm(formula=formula, family=family, data=data)
  }
  #if(varsel!=FALSE)
  if(length(var.rejected)==0){myglm$var.rejected<-"no variables were rejected"
  }else{myglm$var.rejected<-var.rejected#what if no variable selection method is chosen?
  }
  myglm$varsel<- varsel
  if(exists("gn"))myglm$gn<- gn
  myglm$call<- match.call()
  class(myglm)<- "glmsel"
  myglm
}

#' @export
print.glmsel<-function(x, digits=4, ...)#define print method for "linmod class"
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, digits)
  cat("\nNull Deviance:\t  ", format(signif(x$null.deviance, digits)), #signif gets digits significant digits (digits is a devault number)
      "\nResidual Deviance:", format(signif(x$deviance, digits)),
      "\nAIC:\t", format(signif(x$aic, digits)),
      "\nNull df:", x$df.null,
      "\nResid. df: ", x$df.residual)
  if(x$varsel!=FALSE) cat("\nVariable selection method:", x$varsel,
                          "\nVariables rejected:", x$var.rejected)
}

#' Summary of \code{glmsel} objects
#'
#' @aliases summary.glmsel
#' @description This function produces summary of objects of class \code{glmsel} produced by \code{glmsel()} function.
#' @usage \method{summary}{glmsel}(object, dispersion=NULL, ...)
#' @param object an object of class 'glmsel'
#' @param dispersion argument which allows to include dispersion  parameter
#' @param ... arguments to be passed to and from other methods
#' @import stats
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
summary.glmsel<-function(object, dispersion=NULL, ...)#define summary method
{
  #dispersion
  if (is.null(dispersion))
    dispersion <- if (object$family$family %in% c("poisson","binomial"))1
  else if (object$df.residual > 0){
    est.disp <- TRUE
    if (any(object$weights == 0)) warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights * object$residuals^2)[object$weights>0])/object$df.residual
  }
  #dispersion
  #standard errors
  if (object$rank > 0) {
    r <- 1:object$rank
    Qr <- object$qr
    coef.p <- object$coefficients[Qr$pivot[r]]
    covmat.unscaled <- chol2inv(Qr$qr[r, r, drop = FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
    covmat <- dispersion * covmat.unscaled
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
  }
  #standard errors
  tval<-coef(object)/s.err
  TAB<-cbind(Estimate=coef(object),
             Std.Error=s.err,
             t.value=tval,
             p.value=2*pnorm(-abs(tvalue))
  )
  res<-list(call=object$call,
            coefficients=TAB,
            aic=object$aic,
            ndev=object$null.deviance,
            rdev=object$deviance,
            dfnull=object$df.null,
            dfres=object$df.residual,
            varsel=object$varsel,
            var.rejected=object$var.rejected
  )
  class(res)<-"summary.glmsel"
  res
}

#' @export
print.summary.glmsel<-function(x, ...)#adds p.values to the summary method
{#seems it has some restriction on the number of characters printed and doesnt list stars when more vars are present
  #or maybe when intercept is present  it just does not work, problem with has.Pvalue
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.Value=TRUE, has.Pvalue=TRUE)

  cat("\nNull deviance:", x$ndev, "on", x$dfnull, "degrees of freedom",
      "\nResidual deviance", x$rdev, "on", x$dfres, "degrees of freedom",
      "\nAIC:", x$aic
  )
  if(x$varsel!=FALSE) cat("\nVariable selection method:", x$varsel,
                          "\nVariables rejected:", x$var.rejected)
}

#' Plot \code{glmsel} objects
#'
#' @aliases plot.glmsel
#' @description This function produces diagnostic plots of objects of class \code{glmsel} produced by \code{glmsel()} function.
#' @usage \method{plot}{glmsel}(x, ...)
#' @param x an object of class 'glmsel'
#' @param ... arguments to be passed to and from other methods
#' @details \code{plot} will produce residuals versus fitted values and quantile-quantile plots to diagnose the fit of a generalised linear model. If \code{varsel="lasso"} or \code{varsel="enet"} was selected as arguments in \code{lmsel}, \code{plot(object)} will produce and additional plot with lasso or elastic net variable selection paths.
#' @import graphics
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
plot.glmsel<-function(x, ...)
{
  qqnorm(x$residuals)#qq plot for normality of the residuals
  qqline(x$residuals)

  par(mfrow=c(1,1), ask=T)#hit return to see the next plot
  plot(x$fitted.values,x$residuals, xlab="Fitted Values", ylab="Residuals", main="Residuals vs Fitted Values") #resid vs fit plot
  abline(a=0, b=0, lty=2)
  #  plot(x$m.varsel)#
  # plot(x$cv)
  if(x$varsel=="enet")plot(x$gn, label=TRUE, main="Elastic Net variable selection")
  if(x$varsel=="lasso")plot(x$gn, label=TRUE, main="LASSO variable selection")
  par(mfrow=c(1,1), ask=F)

}

#' @name bank
#' @title bank dataset
#' @description the data is related with direct marketing campaigns of a Portuguese banking institution. The marketing campaigns were based on phone calls. Often, more than one contact to the same client was required, in order to access if the product (bank term deposit) would be yes (1) or not (0) subscribed.
#' @docType data
#' @usage data(bank)
#' @format a data frame with 4520 observations on the following 17 variables.
#' @details The data are from Moro (2014) and taken from the UCI ML website \url{https://archive.ics.uci.edu/ml/datasets/Bank+Marketing}.
#' @return \code{bank}	data frame of data with predictor columns \code{age, job, marital, education, default, balance, housing, loan, contact, day, month, duration, campaign, pdays, previous} and \code{poutcome} with response column \code{y} indicating whether the client has subscribed a term deposit.
#' @source S. Moro, P. Cortez and P. Rita. A Data-Driven Approach to Predict the Success of Bank Telemarketing. Decision Support Systems, Elsevier, 62:22-31, June 2014
"bank"

# load plsgenomics library
#library(plsgenomics)
# load leukemia data
#data(leukemia)

#leukemia$X
