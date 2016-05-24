ses.maglm <-
function(data,y,family,scale=TRUE,AIC.restricted=TRUE,par=FALSE, runs=999){

if(scale) data<-as.data.frame(scale(data))

res.obs <- maglm(data,y=y,family,AIC.restricted=AIC.restricted)$importance

# runs<-2#
null.env.list <-list()
# data <- env2 
# before<-proc.time()
for (i in 1:runs){
data.temp <- data
data.temp$temp <-rnorm(nrow(data.temp))
data.temp <- data.temp[order(data.temp$temp),]
null.env.list[[i]] <- data.temp[,-ncol(data.temp)]
}

# y <- "pre.abs"
# family <-"binomial"
# AIC.restricted=F
if (par==FALSE) {res.rand0<-sapply(null.env.list,function(x){mamglm(x,y=y,family,AIC.restricted=AIC.restricted)$importance})}
else res.rand0<-sfSapply(null.env.list,function(x){maglm(x,y=y,family,AIC.restricted=AIC.restricted)$importance})

res.rand<-t(res.rand0) #tranpose (row <-> column)


res.rand.mean <- apply(res.rand,2,mean, na.rm=T)
res.rand.sd <- apply(res.rand,2,sd, na.rm=T)
SES <- (res.obs - res.rand.mean)/res.rand.sd

res.obs.rank <- apply(rbind(res.obs,res.rand),2,rank)[1,]

data.frame(res.obs, res.rand.mean, res.rand.sd, SES, res.obs.rank,runs)

}
