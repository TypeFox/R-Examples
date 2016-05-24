"print.qrssvs"<-function(x, ...){
x.orig<-x
cat("Quantile regression stochastic search \nvariable selection (QR-SSVS) output:\nStart = ", 
	attr(x,"mcpar")[1], "\nEnd = ",
	attr(x,"mcpar")[2], "\nThinning interval = ",
	attr(x,"mcpar")[3], "\n")

attr(x, "mcpar") <- NULL
attr(x, "class") <- NULL
NextMethod("print", ...)
invisible(x.orig)
}

"mptable"<-function(qrssvs){
if (!is(qrssvs, "qrssvs")){
stop("Can only be used on objects of class qrssvs.\n")
}
ssvs.start <- attr(qrssvs, "mcpar")[1]
ssvs.end <- attr(qrssvs, "mcpar")[2]
ssvs.thin <- attr(qrssvs, "mcpar")[3]
nstore <- (ssvs.end-ssvs.start)/ssvs.thin + 1
probs<-apply(qrssvs,2,function(z){length(which(z==1))})/nstore
return(data.frame(Probability=probs))
}

"topmodels"<-function(qrssvs, nmodels=5, abbreviate=FALSE, minlength=3){
if (!is(qrssvs, "qrssvs")){
stop("Can only be used on objects of class qrssvs.\n")
}
ssvs.start <- attr(qrssvs, "mcpar")[1]
ssvs.end <- attr(qrssvs, "mcpar")[2]
ssvs.thin <- attr(qrssvs, "mcpar")[3]
nstore <- (ssvs.end-ssvs.start)/ssvs.thin + 1
xnames <- attr(qrssvs, "xnames")
if (abbreviate){
xnames <- abbreviate(xnames, minlength)
}
model.list<-apply(qrssvs,1,function(z)xnames[which(z==1)])
model.vector<-sapply(model.list, function(z)paste(z, collapse=","))
model.count<-sort(table(model.vector), decreasing=T)/nstore
if (nmodels>length(model.count)){
warning("Number of models requested exceeds total number of models visited.\n")
}
if (rownames(model.count)[1]==""){
rownames(model.count)[1]<-"Null model"
}
return(data.frame(Probability=model.count[1:(min(nmodels, length(model.count)))]))
}

"plot.qrssvs"<-function(x, ...){
probs<-mptable(x)
dotplot(as.matrix(probs),  
  panel=function(x, y, ...){
     panel.abline(v=0.5, lty=3)
     panel.dotplot(x, y, ...)
       },
  origin=0, type=c("p","h"), pch=16, 
  xlim=c(-0.05,1.05),
  scales = list(x = list(at = c(0,0.2,0.4,0.5,0.6,0.8,1))),
  xlab="Marginal inclusion probability", ...)
}

"summary.qrssvs"<-function(object, ...){
covnames <- attr(object, "xnames")
probs<-mptable(object)
median.model<-covnames[probs>=0.5]
results<-probs
attr(results, "median.model")<-median.model
attr(results, "tau")<-attr(object, "tau")
class(results)<-"summary.qrssvs"
return(results)
}

"print.summary.qrssvs"<-function(x, digits=max(3, .Options$digits-3), ...){
attr(x, "class")<-"data.frame"
cat("\nMarginal inclusion probability of each predictor:\n\n")
print(x, digits=digits, ...)
cat("\nFor tau = ", attr(x,"tau"), ", the median probability model \nincludes the following predictors:\n\n",
paste(attr(x, "median.model"), collapse=", "), ".\n\n", sep="")
invisible(x)
}

