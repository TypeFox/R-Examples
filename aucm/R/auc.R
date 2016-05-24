# functions for class auc
# auc has to have fields: formula, coef

ratio <- function(fit) UseMethod("ratio") 
trainauc <- function(fit, ...) UseMethod("trainauc") 

coef.auc=function(object, ...) {
    # for spauc, the field name is coefficient.feature
    tmp=drop(object$coef)
    if (!is.null(names(tmp))) {
        if(names(tmp)[1]=="(Intercept)") {
            tmp=tmp[-1]
        } 
    }
    tmp
}

ratio.auc=function(fit) {
    tmp=coef(fit)
    if (!is.null(names(tmp))) {
        if(names(tmp)[1]=="(Intercept)") {
            tmp=tmp[-1]
        } 
    }
    tmp/tmp[1]
}

ratio.glm=function(fit) {
    tmp=coef(fit)
    if (!is.null(names(tmp))) {
        if(names(tmp)[1]=="(Intercept)") {
            tmp=tmp[-1]
        } 
    }
    tmp/tmp[1]
}


predict.auc=function (object, newdata, case.percentage=NULL, ...) { 
    
    # extract predictor and response
    
    form. = object$formula
    
    y=NULL
    formula.has.response = attr(terms(form.), "response")==1    
    if (formula.has.response) {
        resp.var = as.character(form.)[2]    
        newdata.has.response = resp.var %in% names(newdata)    
        if (newdata.has.response) {
            y=model.frame(form., newdata)[,1]
        } 
    }    
    
    if (formula.has.response & !newdata.has.response) {
        # need to drop response variable from formula
        form. = as.formula(as.character(form.)[-2])
    }
    X.pred = model.matrix(form., newdata)[,-1,drop=FALSE]  
    
    # make predictions
    nonlinear=FALSE
    if (!is.null(object$kernel))
        if (object$kernel!="l") nonlinear=TRUE
    if (nonlinear) {    
    
        if ("rauc" %in% class(object)) {            
            K=getK(X.pred,object$kernel,object$para,X2=object$X) # Note that X.pred should be the first param, not object$X
            Q.pred = getQ(K,n1=object$n.case,n2=object$n.control,call.C=TRUE,do.pred=TRUE)/object$lambda
            res=Q.pred %*% object$alpha.pred
    
        } else if ("dcsauc" %in% class(object) | "srauc" %in% class(object)) {        
            K=getK(X.pred, object$kernel, object$para, X2=object$X) # Note that X.pred should be the first param, not object$X
            res = K %*% object$coefficients
            
        } else stop("predict.auc: something wrong")        
    
    } else {
        res = X.pred %*% ratio(object)
    }        
    res=c(res)
    
    # compute auc
    if (!is.null(y)) {
        attr(res, "auc")=fast.auc(res,y)
    }
        
    # make class prediction
    if (!is.null(case.percentage)) {        
        threshold=quantile(res, 1-case.percentage)
        attr(res,"Class")=as.numeric(res>threshold)
    }
    
    res
}

print.auc=function(x, ...) {
    for (a in names(x)) {
        cat(a,"\n")
        
        if (a %in% c("eta", "linear.combination", "alpha.pred", "y")) {
            print("vector of length "%+%length(x[[a]])%+%" ...", quote=FALSE)
        
        } else if (a=="X" | a=="K") {
            print("matrix of dimension "%+%concatList(dim(x[[a]]), " x ")%+%" ...", quote=FALSE)
        
        } else if (a=="coefficients") {
            if (length(x[[a]])>20) {
                print("coefficients ... vector of length "%+%length(x[[a]])%+%" ...", quote=FALSE)
            } else print(x[[a]], quote=TRUE)
        
        } else if (a=="last.minQuad.fit") {
            print("last minQuad x"%+%" ...", quote=FALSE)            
        
        } else {
            print(x[[a]], quote=TRUE)
        
        }
        cat("\n")
    }
}


summary.auc=function (object, ...) {
    c(trainauc=trainauc(object), minutes=as.double(object$time.elapsed,units="mins"), range=range(object$coefficients))
}


print.minQuad=function(x, quote=TRUE, ...){
    x1=x
    x1$alpha=paste("a numeric vector of size",length(x$alpha))
    class(x1)="list"
    print(x1, quote=quote)
}

trainauc.auc=function(fit, training.data=NULL, ...) {
    if (!is.null(training.data)) {
        y=model.frame(fit$formula, training.data)[,1]
        fast.auc(predict(fit, training.data), y)
    } else if (!is.null(fit$train.auc)) {
        fit$train.auc
    } else NA
} 

trainauc.glm = function (fit, ...) {
    fast.auc(fit$linear.predictors, fit$y)    
}


fast.auc<-function(score, outcome, t0=0, t1=1, reverse.sign.if.nece=TRUE){
    
    if (length(score)!=length(outcome)) stop("score and outcome do not have the same length")
    
   if (all(is.na(score)) | all(is.na(outcome))) {
    
        out=NA
        warning("all score NA or all outcome NA")
    
    } else {
    
        oo<-!is.na(score) & !is.na(outcome)
        if (sum(oo)!=length(score)) warning("some score or outcome are NA")
        if (sum(oo)==0) {
            warning("no non-NA pair")
            return (NA)
        }
        score<-score[oo]; outcome<-outcome[oo]
       
        score1 = score[outcome==1]; n1 = length(score1);
        score2 = score[outcome==0]; n2 = length(score2);
        
        if (t0!=0 | t1!=1) {
            # pAUC
            tmp=quantile(score2, c(1-t1,1-t0))
            score2 = score2[score2>=tmp[1] & score2<=tmp[2]]
        } 
        
        # AUC
        r = rank(c(score1,score2))
        auc = sum(r[1:n1])/n1/n2 - (n1+1)/2/n2
        if (t0==0 & t1==1 & reverse.sign.if.nece) out=max(auc, 1-auc) else out=auc
        
    }
    return(out)
    
}


calAUC<-function(stuff){
#### TPR and FPR must be in decreasing order
    
 TPR<-stuff$sensitivity
 FPR<-1-stuff$specificity
 if ((length(TPR)==1) & (is.na(TPR[1])))
 out<-NA else {
 len<-length(TPR)
 TPRl<-TPR[1:(len-1)]
 TPRu<-TPR[-1]
 diffFPR<-diff(FPR)
 out<-sum((TPRl+TPRu)*diffFPR/(-2))
 }
 return(out)
}
