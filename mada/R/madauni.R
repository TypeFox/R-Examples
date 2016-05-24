
madauni <- function(x, type = "DOR", method = "DSL", suppress = TRUE, ...){

if(suppress){x <- suppressWarnings(madad(x, ...))
             }else{
               x <- madad(x, ...)
             }  
  
TP <- x$data$TP
FP <- x$data$FP
FN <- x$data$FN
TN <- x$data$TN

number.of.pos<-TP+FN
number.of.neg<-FP+TN
nop <- number.of.pos
non <- number.of.neg
total <- nop + non

# from Cochran.Q and inverse variance weights calculate between study variance
naive.tausquared<-function(Q,weights)
{
k<-length(weights)
if(Q<(k-1)){return(0)}
else
return((Q-k+1)/(sum(weights)-(sum(weights^2)/sum(weights))))
}

if(! method %in% c("MH","DSL"))stop("method must be either \"MH\" or \"DSL\"")else
nobs <- x$nobs
theta <-switch(type, "DOR" = x$DOR$DOR, "posLR" = x$posLR$posLR, 
                 "negLR" = x$negLR$negLR)
  
if(method == "MH")
  {
  weights<-switch(type, "DOR" = FP*FN/total, "posLR" = FP*nop/total, 
                  "negLR" = TN*nop/total)
  coef <- log(sum(weights*theta)/sum(weights))
  CQ<-cochran.Q(theta, weights = weights)
  tau.squared <- NULL
  
  P <- sum((nop*non*(TP+FP) - TP*FP*total)/total^2)
  U <- sum(TP*non/total)
  V = sum(FN*nop/total)
  Uprime = sum(FP*non/total)
  Vprime = sum(TN*nop/total)
  R = sum(TP*TN/total)
  S = sum(FP*FN/total)
  E = sum((TP+TN)*TP*TN/(total^2))
  FF = sum((TP+TN)*FN*FP/(total^2))
  G = sum((FP+FN)*TP*TN/(total^2))
  H = sum((FP+FN)*FP*FN/(total^2))
  
  vcov <- switch(type, "DOR" = 0.5*(E/(R^2) + (F+G)/(R*S) + H/(S^2)),
                        "posLR" = P/(U*V), "negLR" = P/(Uprime*Vprime)) 
  }#end of method = "MH"

if(method == "DSL")
  {
  se.lntheta <- switch(type, "DOR" = x$DOR$se.lnDOR, "posLR" = x$posLR$se.lnposLR, 
                       "negLR" = x$negLR$se.lnnegLR)
  lntheta <- log(theta)
  CQ<-cochran.Q(lntheta, 1/se.lntheta^2)
  tau.squared<-naive.tausquared(CQ[1],1/(se.lntheta^2))
  weights<-1/(se.lntheta^2+tau.squared)
  #recalculate CQ based on new weights
  CQ<-cochran.Q(lntheta, weights)  
  coef <- sum(weights*lntheta)/sum(weights)
  vcov <- 1/sum(weights)
  }#end of method == "DSL"

names(coef) <- paste("ln", type, collapse ="", sep ="")

vcov <- matrix(vcov, nrow = 1, ncol = 1)
colnames(vcov) <- paste("ln", type, collapse ="", sep ="")
rownames(vcov) <- paste("ln", type, collapse ="", sep ="")


output <- list(coefficients = coef, vcov = vcov,  tau_sq = tau.squared, weights = weights,
               type = type, method = method, data = x$data, theta = theta, CQ = CQ, nobs = length(theta),
               descr = x, call = match.call())
class(output) <- "madauni"
output
}# end of function madauni
  
print.madauni <- function(x, digits = 3, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if(is.null(x$tau_sq)){
  ans <- exp(x$coefficients)
  names(ans) <- x$type
  print(ans)
      }else{
  ans <- c(exp(x$coefficients),x$tau_sq)
  names(ans) <- c(x$type, "tau^2")
  print(round(ans, digits))
  }  
}

vcov.madauni <- function(object, ...){object$vcov}

summary.madauni <- function(object, level = .95, ...){
x <- object
# Calculate Higgins I^2, only if method is not "MH"  
Higgins.Isq<-function(T,df){return(max(0,100*(T-df)/T))}

if(object$method == "DSL"){
  Isq <- Higgins.Isq(x$CQ[1], x$CQ[3])}else{
  Isq <- NULL
}
  
CIcoef <- rbind(exp(cbind(coef(x), confint(x, level = level))),
                cbind(coef(x), confint(x, level = level)))
rownames(CIcoef) <- c(x$type,paste("ln",x$type, sep ="", collapse = ""))
colnames(CIcoef)[1] <- paste(x$method, "estimate", collapse ="")

Q <- function(tau2){sum(((log(x$theta) - coef(x)[1])^2)/(1/x$weights+tau2))}
CQ <- ifelse(is.null(x$tau_sq), Q(0), Q(x$tau_sq))

## Q-Profile Confidence interval for tau_sq like in Viechtbauer (2007)
if(!is.null(x$tau_sq)){
# browser()
  kappa_up <- qchisq(1-(1-level)/2, x$nobs - 1)
  kappa_low <- qchisq((1-level)/2, x$nobs - 1)
  if(Q(0) < kappa_up){lower <- 0}else{
      lower <-  uniroot(function(x){Q(x)-kappa_up}, lower = 0, upper = 10^4)$root}
  if(Q(0) < kappa_low){upper <- 0}else{
  upper <-  uniroot(function(x){Q(x)-kappa_low}, lower = 0, upper = 10^4)$root}
  CIcoef <- rbind(CIcoef, c(x$tau_sq, lower, upper), sqrt(c(x$tau_sq, lower, upper)))
  rownames(CIcoef)[3:4] <- c("tau^2","tau")
}

output <- list(x=object, Isq = Isq, CIcoef = CIcoef)
class(output) <- "summary.madauni"
output
  
}

print.summary.madauni <- function(x, digits = 3,...){

cat("Call:\n")
print(x$x$call)
cat("\nEstimates:\n")
print(round(x$CIcoef,digits))
cat("\nCochran's Q: ",round(x$x$CQ[1],digits),
    " (",round(x$x$CQ[3])," df, p = ", round(x$x$CQ[2], digits),")", sep = "")
if(!is.null(x$Isq)){cat("\nHiggins' I^2: ",round(x$Isq, digits),"%", sep ="")}
}



