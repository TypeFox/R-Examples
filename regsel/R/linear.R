#devtools::use_package("lars")
#devtools::use_package("glmnet")
#devtools::use_package("MASS")
#devtools::use_package("elasticnet")


linmodEst<-function(x, y) {#linear regression basic function
  qx<-qr(x)#add if terms(myformula)$response==0 then say no response or something
  coef<-solve.qr(qx, y)
  df<-nrow(x)-ncol(x)
  sigma2<-sum((y-as.matrix(x)%*%as.matrix(coef))^2)/df
  vcov<-sigma2*chol2inv(qx$qr)
  colnames(vcov)<-rownames(vcov)<-colnames(x)
  list(coefficients = coef,
       vcov= vcov,
       sigma= sqrt(sigma2),
       df= df,
       qr= qr)
}

#' Linear regression with variable selection
#'
#' \code{lmsel} is used to fit linear models with optionally performing variable selection using stepwise regression, lasso or elastic net methods.
#' @usage lmsel(formula, data=environment(), varsel=FALSE, criterion="AIC",
#' direction="backward", indices=NULL, train=0.3, lambda=1000)
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment containg the variables in the model. If not specified, the variables are taken from the current environment.
#' @param varsel a method of variable selection to be used. The default is \code{"FALSE"}. Available methods include: stepwise regression \code{"step"}, LASSO \code{"lasso"}, elastic net\code{"enet"}.
#' @param criterion when \code{varsel="step", criterion} allows to select a method of calculating statistic for model comparison. The default is  \code{"AIC"}. Less liberal, BIC penalty can be used by typing \code{"BIC"}.
#' @param direction the mode of stepwise search, can be one of \code{"both"}, \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}. If the scope argument is missing the default for \code{direction} is \code{"backward"}.
#' @param indices vector of \code{0} and \code{1} values indicating which observations are to be used as train and test when \code{varsel="lasso" or "enet"}.
#' @param train if \code{indices=NULL}, the function will randomly assign observations as train and test. \code{train} specifies what percentage of data will be used as train observations. Can take values from \code{0.1} to \code{0.9}.
#' @param lambda quadratic penalty parameter for elastic net. The default value is \code{1000}.
#' @return A \code{"lmsel"} object is returned, for which print, plot and summary methods can be used.
#' @examples data(prostate)
#' set.seed(10)
#' lmsel(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45, indices=as.numeric(prostate$train),
#' data=prostate, varsel="lasso")
#'
#' data(concrete)
#' lmsel(CompressiveStrength~., data=concrete, varsel="step", criterion="BIC")
#' @import stats
#' @import elasticnet
#' @import utils
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
lmsel<-function(formula, data=environment(), varsel=FALSE, criterion="AIC", direction="backward", indices=NULL, train=0.3, lambda=1000) UseMethod("lmsel")#define new method #change x into formula??

#' @export
lmsel.default <-function(formula, data=environment(), varsel=FALSE, criterion="AIC", direction="backward", indices=NULL, train=0.3, lambda=1000)#pass lmsel through default method to manipulate and output est
{
  #
  if (!missing(data)) data<- as.data.frame(data)
  x<- model.matrix(formula, data=data)
  y<-model.frame(formula, data=data)[,1]
  if(!is.numeric(y)) stop("response vector is not numeric")
  if(varsel=="enet" & lambda==0) warning("lambda=0 will perform the lasso fit")
  if(train>=1 | train<=0) stop("train must be > 0 and < 1")
  if(lambda<0) stop("lambda cannot be negative")
  if(!all(indices %in% c(0,1))) stop("indices column should only have values 0 and 1")
  if(is.null(indices)==FALSE & all(indices %in% 0)) stop("indices column cannot contain only 0 values")
  if(is.null(indices)==FALSE & all(indices %in% 1)) stop("indices column cannot contain only 1 values")
  #if(missing(indices) & (!missing(varsel))) indices = sample(1:nrow(x), size=train*nrow(x))#will it work if data is unspecified? probably not
  if(is.null(indices) & (varsel!=FALSE)) {indices<-numeric(nrow(x)) #obtain indices for data splitting
  r<-sample(1:nrow(x), size=train*nrow(x))
  indices[r]<-1
  indices[-r]<-0
  }
  var.rejected <- character(0)
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
    m.varsel<- enet(x=x.train, y=y.train,lambda=0)#here x should be train.x #how to choose lambda
    cv<- cv.enet(x=x.train, y=y.train, lambda = 0,s = seq(0,1,0.01), mode="fraction", plot.it=FALSE)#here x should be train.x
    s.enet<-cv$s[which.min(cv$cv)]
    var.all<-predict.enet(m.varsel ,x=x.test,s=s.enet,type="coefficients",mode="fraction")$coefficients#here x should be test.x
    #handle selected factors
    var.all<-var.all[names(var.all)!="(Intercept)"]
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
    x<- model.matrix(formula, data=data)
    y<-model.frame(formula, data=data)[,1]
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
    m.varsel<- enet(x=x.train, y=y.train, lambda=lambda)#here x should be train.x #how to choose lambda
    cv<- cv.enet(x=x.train, y=y.train, lambda = lambda,s = seq(0,1,0.01), mode="fraction", plot.it=FALSE)#here x should be train.x
    s.enet<-cv$s[which.min(cv$cv)]
    var.all<-predict.enet(m.varsel, x=x.test, s=s.enet, type="coefficients",mode="fraction")$coefficients#here x should be test.x
    #handle selected factors
    var.all<-var.all[names(var.all)!="(Intercept)"]
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
    x<- model.matrix(formula, data=data)
    y<-model.frame(formula, data=data)[,1]
  }
  ############
  if(varsel=="step"){
    if(criterion=="BIC") k<-log(nrow(x))
    if(criterion=="AIC") k<-2
    stepest<-step(lm(formula=formula, data=data), k=k, direction=direction, trace=FALSE)
    #extract stepest formula here and should be done
    #x<- model.matrix(formula, data=data)
    #y<-model.frame(formula, data=data)[,1]
    #get formula and data of stepest and plug into linmodEst

    #x<-c("as.factor(season)1", "as.factor(season)2", "hum", "huma", "icos", "as.factor(jakisfactor)g")
    #var.rejected<-c("as.factor(season)", "hum")
    #oddzielnie grep dla factorow i nie factorow, zeby nie bylo usowania hum2 skoro hum jest w var.rehected
    var.rejected<-gsub("^- ", "",(stepest$anova$Step)[-1])

    if(length(var.rejected)!=0)formula<-formula(paste(capture.output(formula),"-",paste(var.rejected, collapse="-")))
    x<-model.matrix(formula, data=data)
  }
  est<-linmodEst(x, y)
  est$fitted.values <-as.vector(x%*%est$coefficients)
  est$residuals<-y-est$fitted.values
  if(exists("m.varsel"))est$m.varsel<-m.varsel
  est$varsel<-varsel
  if(length(var.rejected)==0){est$var.rejected<-"no variables were rejected"
  }else{est$var.rejected<-var.rejected#what if no variable selection method is chosen?
  }
  if(varsel=="lasso" | varsel=="enet") est$cv<- cv#create this object only when lasso or enet are used
  if(varsel=="lasso" | varsel=="enet") est$var.selected<- var
  est$call<-match.call()#create call object  containging function and arguments used
  class(est)<-"lmsel"
  est
}

#' @export
print.lmsel<-function(x, digits=4, ...)#define print method for "lmsel class"
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, digits)
  if(x$varsel!=FALSE) cat("\nVariable selection method:", x$varsel,
                          "\nVariables rejected:", x$var.rejected)
}

#' Summary of \code{lmsel} objects
#'
#' @aliases summary.lmsel
#' @description This function produces summary of objects of class \code{lmsel} produced by \code{lmsel()} function.
#' @usage \method{summary}{lmsel}(object, ...)
#' @param object an object of class 'lmsel'
#' @param ... arguments to be passed to and from other methods
#' @import stats
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
summary.lmsel<-function(object, ...)#define summary method
{
  se<-sqrt(diag(object$vcov))
  tval<-coef(object)/se

  TAB<-cbind(Estimate=coef(object),
             StdErr=se,
             t.value=tval,
             p.value=2*pt(-abs(tval),df=object$df))

  res<-list(call=object$call,
            coefficients=TAB,
            varsel=object$varsel,
            var.rejected=object$var.rejected)
  class(res)<-"summary.lmsel"
  res
}

#' @export
print.summary.lmsel<-function(x, ...)#adds p.values to the summary method
{#seems it has some restriction on the number of characters printed and doesnt list stars when more vars are present
  #or maybe when intercept is present  it just does not work, problem with has.Pvalue
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.Value=TRUE, has.Pvalue=TRUE)

  if(x$varsel!=FALSE) cat("\nVariable selection method:", x$varsel,
                          "\nVariables rejected:", x$var.rejected)

}
#' Plot \code{lmsel} objects
#' @aliases plot.lmsel
#' @description This function produces diagnostic plots of objects of class \code{lmsel} produced by \code{lmsel()} function.
#' @usage \method{plot}{lmsel}(x, ...)
#' @param x an object of class 'lmsel'
#' @param ... arguments to be passed to and from other methods
#' @details \code{plot} will produce residuals versus fitted values and quantile-quantile plots to diagnose the fit of a linear model. If \code{varsel="lasso"} or \code{varsel="enet"} was selected as arguments in \code{lmsel}, \code{plot(object)} will produce and additional plot with lasso or elastic net variable selection paths.
#' @import graphics
#' @author Michal Knut \email{1105406k@student.gla.ac.uk}.
#' @export
plot.lmsel<-function(x, ...)
{
  qqnorm(x$residuals)#qq plot for normality of the residuals
  qqline(x$residuals)

  par(mfrow=c(1,1), ask=T)#hit return to see the next plot
  plot(x$fitted.values,x$residuals, xlab="Fitted Values", ylab="Residuals", main="Residuals vs Fitted Values") #resid vs fit plot
  abline(a=0, b=0, lty=2)
  if(x$varsel=="lasso")plot(x$m.varsel, main="LASSO variable selection")
  if(x$varsel=="enet")plot(x$m.varsel, main="Elastic Net variable selection")
  par(mfrow=c(1,1), ask=F)

}


#' @name prostate
#' @title prostate dataset
#' @description data to examine the correlation between the level of prostate-specific antigen and a number of clinical measures in men who were about to receive a radical prostatectomy.
#' @docType data
#' @usage data(prostate)
#' @format a data frame with 97 observations on the following 10 variables.
#' @details The last column indicates which 67 observations were used as the "training set" and which 30 as the test set, as described on page 48 in the book.
#' @return \code{concrete}	data frame of data with predictor columns \code{lcavol, lweight, age, lbph, svi, lcp, gleason} and \code{pgg45} with response column \code{lpsa} and \code{train} column indicating "training set".
#' @source Stamey, T., Kabalin, J., McNeal, J., Johnstone, I., Freiha, F., Redwine, E. and Yang, N (1989) Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate II. Radical prostatectomy treted patients, Journal of Urology 16: 1076-1083.
"prostate"

#' @name concrete
#' @title concrete dataset
#' @description Yeh (1998) describes a collection of data sets from different sources that can be used for modeling the compressive strength of concrete formulations as a functions of their ingredients and age.
#' @docType data
#' @usage data(concrete)
#' @format a \code{RangedData} instance, 1 row per CpG island.
#' @details The data are from Yeh (1998) and taken from the UCI ML website \url{http://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength}. There are 1030 data points from the UCI website, but the paper states that approximately 1,000 samples were made, but only 727 were analyzed in the source material. It is unclear which samples were excluded.
#' @return \code{concrete}	data frame of data with predictor columns \code{Cement, BlastFurnaceSlag, FlyAsh, Water, Superplasticizer, CoarseAggregate, FineAggregate} and \code{Age} with response column \code{CompressiveStrength}.
#' @source Yeh, I. C. (1998). Modeling of strength of high-performance concrete using artificial neural networks. Cement and Concrete Research, 28(12), 1797-1808. Elsevier.
"concrete"

#var.selected are good, problem with transformation to a formula object, check it manually in liver2 or liver
#add var.selected to generalised or delete it from linear

