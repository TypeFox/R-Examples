## The following is the definition of the thregI function for interval-censored data
"thregI" <- 
function (formula,data) 
{
 ### Read in all arguments
 cl <- match.call()
 indx <- match(c("formula", "data"), names(cl), nomatch=0) 
 if (indx[1] ==0) stop("A formula argument is required")
 mf<- cl[c(1, indx)]
 f <- Formula(formula)
 f1<-formula(f,lhs=1)
 f1<-Formula(f1)
 mf[[1]] <- as.name("model.frame")
 mf$formula <- if(missing(data)) terms(f1) else terms(f1, data=data)  
 mf$formula <- f1
 mf <- eval(mf, parent.frame())
 if (nrow(mf) ==0) stop("No (non-missing) observations")
 Terms <- terms(mf)
 Y <- model.extract(mf, "response")
 if (!inherits(Y, "Surv")) stop("Response must be a survival object")
 type <- attr(Y, "type")
 if (type!='interval') stop(paste("thregI package only support \"", type, "\" survival data", sep=''))
 f2<-formula(f1, lhs = 0)
 if (length(f2[[2]])!=3) stop(paste("Predictors for both lny0 and mu should be specified"))
 x_lny0<-model.matrix(f1, data, rhs=1)
 x_mu<-model.matrix(f1, data, rhs=2)
 left<-Y[,1] 
 right<-Y[,2]
 delta1=matrix(,length(left),1)
 delta2=matrix(,length(left),1)
#right_max=2*max(right[right!=Inf],na.rm=TRUE)
 right_max=10*max(right[right!=Inf])
 for (i in 1 :length(left)){
#if (is.na(right[i])|right[i]=="Inf") {right[i]=right_max
  if (right[i]=="Inf") { right[i]=right_max
                         delta1[i]=0
                         delta2[i]=0}
   else if (left[i]==0) {delta1[i]=1
                         delta2[i]=0}
   else {delta1[i]=0
         delta2[i]=1}
 }

# fix the bug of exact times ................................................................ 
for (i in 1 :length(left))
  {if (right[i]==1 && left[i]>=1) {right[i]=left[i]
                                   left[i]=max(left[i]-0.000001,0)}
  }
 lny0<-function(para_lny0){x_lny0%*%para_lny0}
 mu<-function(para_mu){x_mu%*%para_mu}
 d<-function(para) 
 {
   para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
   para_mu=para[(length(dimnames(x_lny0)[[2]])+1):(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]]))]
   -mu(para_mu)/exp(lny0(para_lny0))
 }
 v<-function(para) 
 {
   para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
   exp(-2*lny0(para_lny0))
 }
 su<-function(para){
   pnorm((1-d(para)*left)/sqrt(v(para)*left))-exp(2*d(para)/v(para))*pnorm(-(1+d(para)*left)/sqrt(v(para)*left))
                   }

 
sv<-function(para){
   pnorm((1-d(para)*right)/sqrt(v(para)*right))-exp(2*d(para)/v(para))*pnorm(-(1+d(para)*right)/sqrt(v(para)*right))
 }
    
 logf<-function(para) {
 -sum(delta1*log(1-sv(para)))-sum(delta2*log(su(para)-sv(para)))-sum((1-delta1-delta2)*log(su(para)))
 }

 p<-rep(0,(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])))
 est<-nlm(logf, p, hessian = TRUE)

 
  names(est$estimate) <-c(paste("lny0:",dimnames(x_lny0)[[2]]),paste("  mu:",dimnames(x_mu)[[2]]))
 loglik = (-1)*est$minimum
 fit<-list(coefficients  = est$estimate,
           var    = solve(est$hessian),
           loglik = loglik,
           AIC    = (-2)*loglik+2*(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])),
           iter   = est$iterations,
           call   = cl,
           mf     = mf,
           lny0   = dimnames(x_lny0)[[2]],
           mu     = dimnames(x_mu)[[2]])
 class(fit) <- 'thregI'
 fit
}