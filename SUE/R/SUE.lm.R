SUE.lm <-
function(formula, data=list(),k, ns, r, constant = 0.25, consistency.check= TRUE){

#define MSE and subsample matrix first
call <- match.call()
   if (missing(data)) 
    data <- environment(formula)
mf=model.frame(formula=formula, data=data)

#compute k,ns,r
if (missing(k)){
para=parameters(nrow(mf), method="appro.k")
k=para$k
r=para$r
ns=para$ns
}

subsample=mf[1:ns,]
for (i in 1:(k-1)){
subsample=rbind(subsample,mf[1:ns,])
}
mse=seq(k)

#fit regression models for k times and generate k MSEs
for (i in 1:k){
choice=sample(1:nrow(mf),ns)
subsample[((i-1)*ns+1):(i*ns),]=mf[choice,]
fit=lm(formula=formula, data=subsample[((i-1)*ns+1):(i*ns),])
mse[i]=anova(fit)$Mean[length(anova(fit)$Mean)]
}

#generate the combined sample
index=order(mse)[1:r]
aa=rbind(subsample[((index[1]-1)*ns+1):(index[1]*ns),],subsample[((index[2]-1)*ns+1):(index[2]*ns),])
Sg=unique(aa)
for (i in 1:(r-2)){
bb=rbind(Sg,subsample[((index[i+2]-1)*ns+1):(index[i+2]*ns),])
Sg=unique(bb)
}

#consistency check
B=seq(k*length(fit$coeff))
dim(B)=c(k,length(fit$coeff))
for (i in 1:k){
fit=lm(formula=formula, data=subsample[((i-1)*ns+1):(i*ns),])
B[i,]=fit$coefficients
}
distance=function(a,b){
t1=abs(a)+abs(b)
t2=1+abs(a)
t1[which(t2>t1)]=t2[which(t2>t1)]
dis=max(abs(a-b)/t1)
dis
}
n=0
d=seq(length(index)-1)
for (i in 1:(length(index)-1)){
d[i]=distance(B[index[1],],B[index[i+1],])
}

#outputs
output=lm(formula=formula, data=Sg)
output$call=call
output$p=list(formula=formula, data=data ,k=k, ns=ns, r=r, constant = 0.25)
output$combined.sample=Sg
output$sample.size=nrow(Sg)
output$mse=sort(mse)[1:length(index)]
output$beta=B[index,]
output$distance=d
output$subsample=subsample
output$MSE=mse
if (any(d > 0.25) == T) output$check= "NO"
else output$check= "YES"
output
}
