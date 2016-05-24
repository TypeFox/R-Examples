waldts <- 
function (object0,object1)
{
if(class(object0)!="LORgee" | class(object1)!="LORgee")
  stop("Both arguments must be objects of 'LORgee' class ")
if(!all(object0$y==object1$y))
  stop("The response variable differs in the two models")
n0 <- length(object0$coefficients)
n1 <- length(object1$coefficients)
if(n0==n1) stop("The two models must be nested")
if(n0<n1) { obj0 <- object0
            obj1 <- object1 } else { 
                   obj0 <- object1
                   obj1 <- object0
                   }
names0 <- names(obj0$coefficients)
names1 <- names(obj1$coefficients)
namestest <- setdiff(names1,names0)
if(length(namestest)==0) stop("The two models must be nested")
if(length(setdiff(names0,names1))!=0) stop("The two models must be nested")
index <- rep(0,length(namestest))
for(i in 1:length(namestest)) index[i] <- which(namestest[i]==names1)
coefftest <- obj1$coefficients[index]
vartest <- obj1$robust.variance[index,index]
waldtest <- t(coefftest)%*%solve(vartest)%*%coefftest 
pvalue <- 1-pchisq(waldtest,length(namestest))
ans <- list(NullModel=obj0$call$formula,AlternativeModel=obj1$call$formula,
            waldstatistic=round(waldtest,4),df=length(namestest),pvalue=round(pvalue,4))
class(ans) <- "waldts"
ans
}