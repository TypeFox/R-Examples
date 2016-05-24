mbes <-
function(formula, data, aux, N=Inf, method="all", level=.95, ...){

### mbes() - model based estimation
# formula = notation for model
# data = data frame with cols id, y, x and N
# aux = known value for mean of secondary information
# method = estimator method (simple,diff,ration,regr,all)
# N = population size N
# level = coverage probability for confidence intervals
# na.action = missing values treatment (now removed), addable for linear model by ...
                                   
### input treatment
if(level<0 | level>1) stop("Wrong input: ", sQuote("level")," has to be probability between 0 and 1.")
if(!is.numeric(N)) stop("Wrong input: ", sQuote("N")," is not a number or ", sQuote("Inf"),".")     
if(!(method=="simple" || method=="diff" || method=="ratio" || method=="regr" || method=="all")){
  stop("\n Wrong input: Only methods ", sQuote("simple"),", ", sQuote("diff"),", ", sQuote("ratio"),", ", sQuote("regr")," or ", sQuote("all")," are allowed.\n")
}
# request formula
if(missing(formula)) stop("Wrong input: Missing formula.")
else{
  varNames <- attr(terms(formula),'term.labels')
  p <- length(varNames)
}
# request data
if(missing(data)){
  for(i in 1:p){
   if(!exists(varNames[i])){  
     stop("Wrong input: Missing data or wrong input of data")
   }
  }
  data <- model.frame(formula, na.action=na.omit)
}
# data and formula available
else{   
  if(!(is.matrix(data)||is.data.frame(data))) stop("Wrong input: data must be of type matrix or data.frame")
  data <- as.data.frame(data)
}
# match data and formula
varPos <- match(varNames, colnames(data))
if(any(is.null(varPos))) stop("Wrong input: formula and data mismatch")
mf <- model.frame(formula, data=data, na.action=na.omit)
n <- nrow(mf)

# check auxiliary information
if(missing(aux)) stop("Error: If there is no true mean of auxiliary variable ", sQuote("mbes()")," does not suit. Use function ", sQuote("Smean()"),".")
else{
  if(length(aux) > 1 && method!="regr") stop("Wrong input: More than one auxiliary variable just allowed for method=", sQuote("regr"),".")
  if(length(aux) != p) stop("Wrong input: Length of ", sQuote("aux")," has not a value for every auxiliary variable.")
}
### return argument
  ret <- list()
  ret$call <- list(formula=formula,data=data,aux=aux,N=N, method=method, level=level)
  x.mean <- ifelse(ncol(mf)>2,colMeans(mf[,-1], na.rm=TRUE), mean(mf[,-1],na.rm=TRUE))
  ret$info <- list(N=N,n=n,p=p,aux=aux,x.mean=x.mean)
  ret$simple <- list()
  ret$diff <- list()
  ret$ratio <- list()
  ret$regr <- list()

### simple sampling
if (method=="simple"|method=="all"){
  est <- mean(mf[,1], na.rm=TRUE)
  if (! is.infinite(N))
  var.est <- var(mf[,1])/n * (N-n)/N
  else var.est <- var(mf[,1])/n
  # confidence
  kr <- est + qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  kl <- est - qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  if (abs(kr)<1) kr<- signif(kr, 4) else kr<- round(kr, 4)
  if (abs(kl)<1) kl<- signif(kl, 4) else kl<- round(kl, 4)
  # calculation of se
  if(var.est>0) se <- sqrt(var.est)
  else{
    se <- NA
    kl <- NA
    kr <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals cannot be calculated.")
  }  
  # return argument
  ret$simple <- list(mean=est,se=se, ci=c(kl,kr))
}
### diff est
if (method=="diff"|method=="all"){
  est <- mean(mf[,1],na.rm=TRUE) + aux - mean(mf[,2],na.rm=TRUE)
  d <- mf[,1]-mf[,2]
  if (! is.infinite(N))
  var.est <- var(d)/n * (1- n/N)
  else var.est <- var(d)/n
  # confidence
  kr <- est + qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  kl <- est - qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  if (abs(kr)<1) kr<- signif(kr, 4) else kr<- round(kr, 4)
  if (abs(kl)<1) kl<- signif(kl, 4) else kl<- round(kl, 4)
  # calculation of se
  if(var.est>0) se <- sqrt(var.est)
  else{
    se <- NA
    kl <- NA
    kr <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals cannot be calculated.")
  }  
  # return argument
  ret$diff <- list(mean=est,se=se, ci=c(kl,kr))
}
### ratio est
if (method=="ratio"|method=="all") {
  b <- mean(mf[,1], na.rm=TRUE)/mean(mf[,2], na.rm=TRUE)
  est <- mean(mf[,1], na.rm=TRUE) + b*(aux-mean(mf[,2], na.rm=TRUE))
  if (!is.infinite(N))
  var.est <- (N-n)/N * 1/(n*(n-1))*(sum((mf[,1]-b*mf[,2])^2) )
  else var.est <- 1/(n*(n-1))*(sum((mf[,1]-b*mf[,2])^2) )
  # confidence
  kr <- est + qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  kl <- est - qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  if (abs(kr)<1) kr<- signif(kr, 4) else kr<- round(kr, 4)
  if (abs(kl)<1) kl<- signif(kl, 4) else kl<- round(kl, 4)
  # calculation of se
  if(var.est>0) se <- sqrt(var.est)
  else{
    se <- NA
    kl <- NA
    kr <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals cannot be calculated.")
  }  
  # return argument
  ret$ratio <- list(mean=est,se=se, ci=c(kl,kr))
}
### regression est
if (method=="regr"|method=="all") {
  b <- lm(formula, data, ...)
  if(any(is.na(b$coefficients)) ){
    stop("\nLinear Regression Model:",print(b),"\n", sQuote("NA")," coefficients detected. Please remove entry nr.", which(is.na(b$coefficients)))
  }
  if(attr(terms(formula),'intercept')==1)  est <- b$coefficients[1] + sum(b$coefficients[-1] * aux )
  else  est <- sum(b$coefficients * aux )
# calculate var.estimate
  E2 <- sum(b$residuals^2)
  df <- b$df.residual
  if (N<Inf){
    var.est <- (N-n)/N * 1/(n*df)*E2
  }
  else{
    var.est <- 1/(n*df)*E2
  }
  # confidence boundaries
  kr <- est + qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  kl <- est - qnorm(level+(1-level)/2,0,1)*sqrt(var.est)
  if (abs(kr)<1) kr<- signif(kr, 4) else kr<- round(kr, 4)
  if (abs(kl)<1) kl<- signif(kl, 4) else kl<- round(kl, 4)
  # calculation of se
  if(var.est>0) se <- sqrt(var.est)
  else{
    se <- NA
    kl <- NA
    kr <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals cannot be calculated.")
  }  
  # return argument
  ret$regr <- list(mean=est,se=se, ci=c(kl,kr), model=b)
}
structure(ret,class="mbes")
}
