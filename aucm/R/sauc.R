# logistic approximation of SAUC
logistic.f=function(eta,h,loss=TRUE){
    tmp=eta/h
    if (loss) tmp=-tmp
    drop( expit(tmp) )
}


# Normal CDF approximation of SAUC
phi.f=function(eta,h,loss=TRUE){
    out=1-pnorm(eta/h,0,1)
    if (loss) drop(out) else drop(1-out)
}


# from Ying and Shuxin
# not up to usual standard, e.g. no convergence attribute in the return object
# constrain.method="L2"; h.method="Lin"; start.method="rlogit"; opt.method="YH"; upper=NULL  #default
sauc.phi<-function(formula,dat,constrain.method="L2",h.method="Lin",start.method="rlogit", opt.method="Lin", upper=NULL, verbose=FALSE, ret.vcov=FALSE, truth=NULL, beta.init=NULL) {
    
    tmp=model.frame(formula, dat)
    Y=tmp[,1]
    n1=sum(Y==1)
    n2=sum(Y==0)
    n=n1+n2
    predictors=model.matrix(formula, dat)[,-1] # assume there is an intercept in the formula
    pp=ncol(predictors)
    X1=model.matrix(formula, tmp[tmp[,1]==1,])[,-1,drop=FALSE]
    X2=model.matrix(formula, tmp[tmp[,1]==0,])[,-1,drop=FALSE]
    
    if (is.null(beta.init)) {
        if (start.method=="rlogit") {
            fit.rlogit=rlogit(formula, dat)
            if (fit.rlogit$convergence) {
                beta.init=fit.rlogit$coef[-1]
            } else {
                beta.init=rep(1, ncol(tmp)-1)
            }        
        } else if (start.method=="1") {
            beta.init=rep(1, ncol(tmp)-1)
        } else if (start.method=="0") {
            beta.init=rep(0, ncol(tmp)-1)
        } else stop("start.method not supported: "%+%start.method)
    }
    
    Xint<-calXdiff(predictors,Y,d=length(beta.init))
    x_diff<-Xint[[1]];
    x_diff_mul<-Xint[[2]]
    
    if (h.method=="Lin") {
        if (constrain.method=="L2") {
            h=n^(-1/3)*sd(drop(x_diff%*%normsq(beta.init)))
        } else if (constrain.method=="beta1"){
             h=n^(-1/3)*sd(drop(x_diff%*%normanc(beta.init)))
        } else stop("constrain.method not supported: "%+%constrain.method)
    } else if (h.method=="MH") {
        # Ma and Huang
        if (constrain.method=="L2") {
            h=min(1/sqrt(sum(Y==1)),1/sqrt(sum(Y==0)),quantile(abs(x_diff%*%normsq(beta.init)),0.05)/5)        
        } else if (constrain.method=="beta1"){
            h=min(1/sqrt(sum(Y==1)),1/sqrt(sum(Y==0)),quantile(abs(x_diff%*%normanc(beta.init)),0.05)/5)
        } else stop("constrain.method not supported: "%+%constrain.method)
    } else if (h.method=="Vexler") {
        h=(n1*n2)^(-0.1)
    } else if (is.numeric(h.method)) {
        h=h.method
    } else stop("h.method not supported: "%+%h.method)
    
    
    
    if (opt.method=="truth") {
        out = truth
    
    } else if (opt.method=="YH") {
        if (ncol(x_diff)==2) {
            if (is.null(upper)) stop("upper has to be supplied when opt.method is YH")
            if (constrain.method=="L2") {
                tmp=cal.SGD(beta.init,x_diff, upper=upper, h, verbose)
                out=tmp[-length(tmp)]
            } else if (constrain.method=="beta1") {
                tmp=calAnchor.SGD(beta.init,x_diff, upper, h, verbose)
                out=tmp[-length(tmp)]
            } else stop("constrain.method not supported: "%+%constrain.method)
        } else stop ("opt.method YH only implemented for two covariates")
        
    } else {
        if (constrain.method=="L2") {        
            beta_opt=try.error(cal(normsq(beta.init),x_diff,x_diff_mul,rep(1,nrow(x_diff)),h))
            #print(beta_opt)
            if(!inherits(beta_opt,'try-error')) {        
                auc.end<-auc_aem(x_diff,normsq(beta_opt[1:pp]),rep(1,nrow(x_diff)),h)
                auc.init<-auc_aem(x_diff,normsq(beta.init),rep(1,nrow(x_diff)),h)
                
                #print(c(auc.init,auc.end))
                out=beta_opt[1:pp]
                #   if (!is.na(auc.init) & !is.na(auc.end) & auc.end<auc.init) out<-rep(NA,2)
            }  else out=rep(NA,2)
            
        } else if (constrain.method=="beta1") {
        
            beta_opt=try.error(calAnchor(normanc(beta.init),x_diff,x_diff_mul,rep(1,nrow(x_diff)),h))
            if(!inherits(beta_opt,'try-error')) {        
                auc.end<-auc_aem(x_diff,normanc(beta_opt[1:pp]),rep(1,nrow(x_diff)),h)
                auc.init<-auc_aem(x_diff,normanc(beta.init),rep(1,nrow(x_diff)),h)
                #print(c(auc.init,auc.end))
                out<-beta_opt[1:pp]
                #  if (!is.na(auc.init) & !is.na(auc.end) & auc.end<auc.init) out<-rep(NA,2)
            }  else out<-rep(NA,2)
            
        } else stop("constrain.method not supported: "%+%constrain.method)
    
    }    
    
    res=list()
    class(res)=c("auc","sauc",class(res))
    res$coefficients=out
    res$formula=formula
    res$kernel="l"
    
    res$train.auc = fast.auc( c(X1 %*% out, X2 %*% out), c(rep(1, n1), rep(0, n2)) )
    
    res$h=h
    if (ret.vcov) res$vcov = get.vcov.sauc ("normal", beta.hat=res$coefficients/res$coefficients[1], h=h/res$coefficients[1], X1, X2)
    
    return(res)
}

get.vcov.sauc=function(sigm.f, beta.hat, h, X1, X2) {
    
    if (!sigm.f %in% c("logistic","normal")) stop("this sigmoid approximation not supported")
    coef.=c(beta.hat)
    d=ncol(X1)
    n1=nrow(X1)
    n2=nrow(X2)
    n=n1+n2
    
    tmp = sapply(1:n1, simplify="array", function(i) {
        ta = switch (sigm.f,
            logistic = {
                aux = expit ( - coef. %*% (X1[i,] - t(X2)) / h); # to use expit, we need a minus sign in there
                aux * (1-aux)
            }, 
            normal = dnorm (coef. %*% (X1[i,] - t(X2)) / h)  
        )
        ta1=rowMeans(rep.matrix(ta,d-1) * (X1[i,] - t(X2))[-1,])
        ta1=as.vector(ta1)
        outer(ta1, ta1)
    })
    tmp2 = sapply(1:n2, simplify="array", function(j) {
        ta = switch (sigm.f,
            logistic = {
                aux  =expit ( - coef. %*% (t(X1) - X2[j,]) / h);
                aux * (1-aux)
            },
            normal = dnorm (coef. %*% (t(X1) - X2[j,]) / h)
        ) 
        ta1 = rowMeans(rep.matrix(ta,d-1) * (t(X1) - X2[j,])[-1,])
        ta1=as.vector(ta1)
        outer(ta1, ta1)
    })
    if (is.null(dim(tmp))) {
        T1=mean(tmp)*(n1/n)*(n2/n)^2 + mean(tmp2)*(n1/n)^2*(n2/n)
    } else {
        T1=apply(tmp, 1:2, mean)*(n1/n)*(n2/n)^2 + apply(tmp2, 1:2, mean)*(n1/n)^2*(n2/n)
    }
    T1 = T1 / h^2
    
    X.diff = get.X.diff(X1, X2)
    tmp2 = switch (sigm.f,
        logistic = {
            tmp = expit( - X.diff %*% coef./h )
            c(tmp*(1-tmp)*(1-2*tmp))
        },
        normal = {
            tmp=X.diff %*% coef./h
            c(tmp*dnorm(tmp))
        }
    )
    T2 = crossprod(X.diff[,-1] * tmp2, X.diff[,-1]) /n^2 /h^2
    T2.inv=solve(T2)   
    
    T2.inv %*% T1 %*% T2.inv / n # divide by n to get variance of beta.hat, otherwise it is the asymptotic variance
    
    # for debugging purpose
    #diag(c(T1, T2, sqrt(T2.inv %*% T1 %*% T2.inv / n))) # diag is here to satisfy downstream code
    
}


vcov.sauc=function(object, ...) {
    object$vcov
}

try.error <- function(expr, silent=TRUE) {

# make it to error
    op <- options("warn")
    on.exit(options(op))
    options(warn = 1)

# catch using the try() function

    try(expr, silent)
    }



normsq<-function(x) {
    return(x/sqrt(sum(x^2)))
}

normanc<-function(x) {
    return(x/x[1])
}

#Calculate empirical AUC
auc_em=function (x,beta,w) {
    ind<-(x%*%beta>0)+(x%*%beta==0)/2
    auc=1/(sum(w))*sum(w*ind)
        result=max(auc,1-auc)
               return(result)
}

#Calculate approximated empirical AUC
auc_aem=function (x,beta,w,h) {
    auc=1/(sum(w))*sum(w*pnorm(x%*%beta/h))
        result=max(auc,1-auc)
               return(result)
}


dauc_aem=function (x,beta,w,h) {
    auc=1/(sum(w))*sum(w*dnorm(x%*%beta/h)*x[,2]/h)
        return(auc)
}




d2auc_aem=function (x,beta,w,h) {
    auc=1/(sum(w))*sum(w*dnorm(x%*%beta/h)*x[,2]/h*(-1)*(x%*%beta/h)*x[,2]/h)
        return(auc)
}

###First derivative of function S
s_der_1_w=function (beta,x,w,h) {
    value=x*as.vector(dnorm(x%*%beta/h)*w/h/dim(x)[1])
          result=apply(value,2,sum)
                 return(result)
}

###Second derivative of function S
s_der_2_w=function (beta,x,x_mul,w,h) {
    value=-x_mul*as.vector(dnorm(x%*%beta/h)*(x%*%beta)/h^3*w/dim(x)[1])
          result=matrix(apply(value,2,sum),length(beta),length(beta),byrow=TRUE)
                 return(result)
}

#Calculate Optimal Beta Starting with pair difference ...
cal=function(beta,x,x_mul,w, h) {
    beta_hat<-rep(NA,length(beta))
    beta.old=beta
             iter=0
    repeat {
        iter=iter+1
        s1=s_der_1_w(beta.old,x,w,h)
        s2=s_der_2_w(beta.old,x,x_mul,w,h)
        beta.new=beta.old-solve(s2)%*%s1
        beta.new=beta.new/sqrt(sum(beta.new^2))
        if (sum(abs(beta.new-beta.old))<5e-4 | iter>200) break
        beta.old=beta.new
    }
    beta_hat=beta.new
    return(c(beta_hat,iter))
}


calAnchor<-function(beta,x,x_mul,w, h) {
#### fix the first beta to be 1  ####
    beta_hat<-rep(NA,length(beta))
    beta.old=beta
             iter=0
    repeat {
        iter=iter+1
        s1=s_der_1_w(beta.old,x,w,h)
        s2=s_der_2_w(beta.old,x,x_mul,w,h)
        beta.new=beta.old[-1]-solve(s2[-1,-1])%*%s1[-1]
        beta.new=c(1,beta.new)
        if (sum(abs(beta.new-beta.old))<5e-4 | iter>200) break
            beta.old=beta.new
        }
    if (iter<=200) beta_hat=beta.new
                                return(c(beta_hat,iter))
    }



### generate pair difference ...
calXdiff<-function(X,Y,d) {
    x_nd<-X[Y==0,]
    x_d<-X[Y==1,]

    n1=sum(Y==1);
    n0=sum(Y==0)
       x_ndr= apply(x_nd, 2, rep, n1)
              x_dr<-apply(x_d,2,rep,each=n0)

              x_diff=x_dr-x_ndr


                     x_diff_mul=matrix(0,dim(x_diff),d*d)
for (p in 1:d) {
for (q in 1:d) {
            x_diff_mul[,(p-1)*d+q]=as.matrix(x_diff[,p]*x_diff[,q])
        }
    }
                           out=list(x_diff,x_diff_mul)
                               return(out)
}


# Ying's gradient based optimization method
cal.SGD=function(beta,x, upper, h, verbose=FALSE){
    beta_hat<-rep(NA,length(beta))
    beta.old=normsq(beta)
    iter=0
    repeat {
     iter=iter+1
     s1=s_der_1_w(beta.old,x,rep(1,nrow(x)),h)
     f<-function(lambda){
         mean(pnorm(x%*%(normsq(beta.old+lambda*s1))/h))
     }
    
    op=optimize(f,c(0,upper),maximum=TRUE)
    lambda=op$maximum     
     beta.new=beta.old+lambda*s1
     beta.new=normsq(beta.new)
     if (sum(abs(beta.new-beta.old))<1e-5 | iter>200) break
     beta.old=beta.new
    if (verbose) print(c(beta.new,op$objective))
    }
    
    if (iter<=200) beta_hat=beta.new
    return(c(beta_hat,iter))
}
#Calculate Optimal Beta Starting with pair difference ... 
calAnchor.SGD=function(beta,x, upper, h, verbose=FALSE){
    beta_hat<-rep(NA,length(beta))
    beta.old=normanc(beta)
    iter=0
    repeat {
     iter=iter+1
     s1=s_der_1_w(beta.old,x,rep(1,nrow(x)),h)
     f<-function(lambda){
         mean(pnorm(x%*%(normanc(beta.old+lambda*c(0,s1[-1])))/h))
     }
    
    op=optimize(f,c(0,upper),maximum=TRUE)
    lambda=op$maximum     
     beta.new=beta.old+lambda*c(0,s1[-1])
     beta.new=normanc(beta.new)
     if (sum(abs(beta.new-beta.old))<1e-5 | iter>200) break
     beta.old=beta.new
     if (verbose) print(c(beta.new,op$objective,op$maximum))
    }
    if (iter<=200) beta_hat=beta.new
    return(c(beta_hat,iter))
    
}
