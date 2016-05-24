
setClass("bild", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric", residuals="numeric", 
		s.residuals="numeric",ind.probability="numeric", prob.matrix="matrix", Fitted="numeric", 
		Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame", y.av="numeric", f.value="factor",call="language"))


setMethod(f="anova", signature(object = "bild"), 
function(object,...,test=TRUE,correct=FALSE)
{

dots <- list(...)

object <- list(object, ...)

if(length(object)<2)
stop("single argument anova not implemented")

x<-object[[1]]
y<-object[[2]]

data1<-x@n.cases
data2<-y@n.cases

if(data1!=data2)
stop("all models must be fit to the same data object")

m1<-x@call$formula
m2<-y@call$formula
x3<-x@call$dependence
y3<-y@call$dependence

if(m1==m2&&x3==y3 )
stop("models are identical")


if(length(x@coefficients)>length(y@coefficients)){

x1<-x@log.likelihood
x2<-x@aic 

 	
y1<-y@log.likelihood
y2<-y@aic 

n1<-length(x@coefficients)
n2<-length(y@coefficients)

n11<-length(x@Time)
n22<-length(y@Time)

bic1<--2*x@log.likelihood+length(x@coefficients)*log(n11)
bic2<--2*y@log.likelihood+length(y@coefficients)*log(n22)

X3<-c(bic1,bic2)

df<-length(x@coefficients)-length(y@coefficients)

}
else {

x1<-y@log.likelihood
x2<-y@aic 
 	
y1<-x@log.likelihood
y2<-x@aic 

data1<-y@n.cases
data2<-x@n.cases


n11<-length(y@Time)
n22<-length(x@Time)


df<-length(y@coefficients)-length(x@coefficients)

bic1<--2*y@log.likelihood+length(y@coefficients)*log(n11)
bic2<--2*x@log.likelihood+length(x@coefficients)*log(n22)


X3<-c(bic2,bic1)

}

X1<-c("Model1  ","Model2  ")
X2<-c(x@aic,y@aic)

X31<-c(x@log.likelihood,y@log.likelihood)

if(test==TRUE){

## Test. L. Ratio

  teste<-2*(x1-y1)# x1 mais geral

if(correct==FALSE){
  p<-1-pchisq(teste,df)}
else{

p<-0.5*(1-pchisq(teste,df))

}


X4<-c(" ",round(teste,3))
X5<-c(" ",round(df,0))
X6<-c(" ",formatC(p))

tabela1<-data.frame(X1,X2,X3,X31,X4,X5,X6)
names(tabela1)<-c(" ","AIC","BIC","logLik"," Deviance","df", "p-value")

cat("\nData:  ")
print(x@call$data)

cat("\nModel1:  ")
print(x@call$formula)
a<-x@call$dependence
cat("dependence =",a,"\n")


cat("Model2:  ")
print(y@call$formula)
b<-y@call$dependence
cat("dependence =",b,"\n")
print(tabela1,row.names=FALSE)}

else{
tabela2<-data.frame(X1,X2,X3,X31)

names(tabela2)<-c(" ","AIC","BIC","logLik")

cat("\nData:  ")
print(x@call$data)

cat("\nModel1:  ")
print(x@call$formula)
a<-x@call$dependence
cat("dependence =",a,"\n")


cat("Model2:  ")
print(y@call$formula)
b<-y@call$dependence
cat("dependence =",b,"\n")
print(tabela2,row.names=FALSE)}


}
)

