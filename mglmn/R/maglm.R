maglm <-
function(data,y,family,scale=TRUE,AIC.restricted=FALSE){
if(scale) data<-as.data.frame(scale(data))
my.vars <- colnames(data)
n.vars <- length(my.vars)

vars.list<-list()
for (i in 1:n.vars){
vars.list[[i]] <-combn(my.vars,i)
}

#numbers of combination for each sample size
temp<-sapply(vars.list,ncol)

#AIC for each model
model.aic <-NULL
log.L<-NULL
vars<-list()
k <- 0
vars2<-matrix(numeric(n.vars*sum(temp)),nrow=sum(temp),ncol=n.vars)

#number of paremeter = i+1 (one means intercept)
#sample size = nrow(data)
#AICc = -2*logLike + 2*n.para*n.sample/(n.sample-n.para-1)
for (i in 1:n.vars){
for (j in 1:temp[i]){
k <- k + 1
vars.temp <- vars.list[[i]][,j]
vars[[k]] <- vars.temp
f.str <- make.formula(y, vars.temp)
# if (family=="gaussian") fit.temp <- manylm(f.str,data=data)
# else 
fit.temp <- glm(f.str,data=data, family=family)

log.L.temp <- logLik(fit.temp)
log.L <- c(log.L, sum(log.L.temp))

if (AIC.restricted){
#aic is corrected for finite sampe size
aic.restricted <- sum(-2*log.L.temp +2*nrow(data)*(i+1)/(nrow(data)-(i+1)-1))
# aic.restricted <- sum(-2*log.L.temp +2*nrow(data)*(i+1)/(nrow(data)-(i+1)-1))
model.aic <- c(model.aic, aic.restricted)
} 
else{
# aic.unrestricted <- sum(-2*log.L.temp +2*(i+1))
aic.unrestricted <- sum(-2*log.L.temp +2*(i+1))

model.aic <- c(model.aic, aic.unrestricted)
} 

}
}

min.aic <- min(model.aic)
delta.aic <- model.aic - min.aic
wAIC <- exp(-delta.aic/2)/sum(exp(-delta.aic/2))
res <- data.frame(AIC=model.aic, log.L=log.L,delta.aic, wAIC,n.vars=rep(1:n.vars,temp))

##counting vars
#vars2 is matrix filled with 0 (row:sites,col:parameters)
vars2<-matrix(numeric(n.vars*sum(temp)),nrow=sum(temp),ncol=n.vars)
colnames(vars2) <- my.vars
n.size<-rep(1:n.vars,temp) # number of paremters for each row
for (i in 1:nrow(vars2)){ # each row (sample)
for (j in 1:n.vars){ # each column (paramter)
for (k in 1:n.size[i]){ # upto number of paramters
if (colnames(vars2)[j]==vars[[i]][k]) vars2[i,j] <- 1
}
}
}

# for exmaple, i=100
# > vars[[100]]
# [1] "LOG..x.1..Exotic.plant.rich"          
# [2] "P.mega.Is.p.a"                        
# [3] "Tot.rain.during.sampling.log..x.0.01."
# these three paramters were used in the 100th analysis
#
# if jth colnames of vars2, which is a parameter name, is identical to vars[[i]][k], which is a paramter name used in the analysis, vars2[i,j] is replaced by 1

res <-cbind(res, vars2)
res <- res[order(res$AIC),]

rownames(res) <- NULL

#calculating weighted result of explanable variables
res.temp <- res[,-1:-5]
res2<-apply(apply(res.temp, 2, function(x)res$wAIC*x),2,sum)
list(res.table=res,importance=res2,family=family)
}
