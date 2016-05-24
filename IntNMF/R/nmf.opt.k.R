nmf.opt.k <-
function(dat=dat,n.runs=30,n.fold=5,k.range=2:8,result=TRUE,make.plot=TRUE,progress=TRUE,
st.count=10,maxiter=100)
{
if (!is.list(dat)) dat <- list(dat) 
n.dat <- length(dat) 
for (i in 1:n.dat)  assign(paste("d",i,sep=""),eval(parse(text=paste("dat[[",i,"]]",sep=""))))
for (i in 1:n.dat){
if (!all(eval(parse(text=paste("d",i,sep="")))>=0))
stop(paste("All values must be positive. There are -ve entries in dat",i,sep=""))}
set.seed(12345)                         
n.sample <- nrow(dat[[1]])

CPI <- matrix(NA,length(k.range),n.runs)   # Adj Rand index between predicted and computed test sample cluster membership
dimnames(CPI) <- list(paste("k",k.range,sep=""),paste("run",1:n.runs,sep=""))

count <- 0
for (i in 1:n.runs){
for (k in k.range){
## Cluster Prediction Index : Agreement between the predicted and computed cluster membership for test data
## Use k-fold cross validation
R.ind <- NULL
random.sample <- sample(seq(n.sample),n.sample)
for (j in 1:n.fold){
test.sample <- random.sample[(round((j-1)*n.sample/n.fold) +1): round(j*n.sample/n.fold)]  
train.sample <- setdiff(random.sample,test.sample) 

d.train <- lapply(dat, function(x) x[train.sample, ,drop=FALSE]) 
d.test <- lapply(dat, function(x) x[test.sample, ,drop=FALSE])
# Estimate Hi's from the training data, compute W using the test data and 
# find cluster membership for the test data  
fit.train <- nmf.mnnals(dat=d.train,k=k,maxiter=maxiter,st.count=st.count,n.ini=1,ini.nndsvd=FALSE,seed=FALSE)
for (m in 1:n.dat){ 
if(n.dat > 1){
assign(paste("H.train",m,sep=""),eval(parse(text=paste("fit.train$H[[",m,"]]",sep=""))))
} else {
assign("H.train",eval(parse(text="fit.train$H")))}
}
XHt <- 0
for (m in 1:n.dat){
if(n.dat > 1){
XHt <- XHt + d.test[[m]]%*%t(eval(parse(text=paste("H.train",m,sep=""))))
} else {
XHt <- XHt + d.test[[m]]%*%t(eval(parse(text="H.train")))}
}
HHt <- 0
for (m in 1:n.dat){ 
if(n.dat > 1){
HHt <- HHt + eval(parse(text=paste("H.train",m,sep="")))%*%t(eval(parse(text=paste("H.train",m,sep="")))) 
} else {
HHt <- HHt + eval(parse(text="H.train")) %*% t(eval(parse(text="H.train")))}
}
W.predict <- XHt %*% ginv(HHt)
W.predict <- W.predict + abs(min(W.predict))
predicted.cluster.mem <- apply(W.predict,1,which.max)

# Apply integrative NMF to test data to compute cluster membership
fit.test <- nmf.mnnals(dat=d.test,k=k,maxiter=maxiter,st.count=st.count,n.ini=1,ini.nndsvd=FALSE,seed=FALSE)
computed.cluster.mem <- fit.test$clusters 

# Adjusted rand index between the predicted and computed cluster membership
R.ind <- c(R.ind,adjustedRandIndex(predicted.cluster.mem,computed.cluster.mem))
# Display progress %
count <- count + 1
if(progress & round(count/(n.runs*length(k.range)*n.fold)*100, 0) %in% seq(5,100,by=5)){ 
message(paste(round(count/(n.runs*length(k.range)*n.fold)*100, 0),"% complete",sep=""))
flush.console()}
}
CPI[k-1,i] <- mean(R.ind)
}
}
# Plots of optimality criteria vs k
if(make.plot){
dev.new(width=4, height=5)
plot(k.range,CPI[,1],ylim=c(min(CPI),max(CPI)),pch=20,main="",xlab="k",ylab="CPI")
for (m in 2:n.runs) points(k.range,CPI[,m],pch=20)
lines(k.range,apply(CPI,1,mean),col="red",lwd=2)
mtext("Optimum k", outer = TRUE, cex = 1, line=-2)}
if(result) return(CPI) 
}
