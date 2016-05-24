
MRsurvival <- function(time.ref, status.ref, cov.rel, data.rel, cox.ref, cov.ref,  init, B)
{

coef.simul <- data.frame(rmvnorm(B, mean=cox.ref$coefficient, sigma=cox.ref$var))
rel.data <- data.rel
n <- dim(rel.data)[1]
rel.data$ident <- 1:n
res.theta <- matrix(-99, nrow=B, ncol=length(cov.rel))
dimnames(res.theta)[[2]] <- cov.rel

for (b in 1:B)
{
num<-sample(rel.data$ident, size=n, replace = TRUE)
data.boot<-rel.data[num,]
data.boot<-data.boot[order(data.boot[,time.ref]),]

data.boot$nb<-1
for (i in unique(data.boot[data.boot[,status.ref], time.ref])) {
data.boot$nb[data.boot[,time.ref]==i & data.boot[,status.ref]==1] <- 1:sum(data.boot[,time.ref]==i & data.boot[,status.ref]==1) }
data.boot$Tps.Evt.Cor<-data.boot[,time.ref]+data.boot$nb/1000
data.boot<-data.boot[order(data.boot$Tps.Evt.Cor, 1-data.boot[,status.ref]),]
data.boot$indic<-1:n

covarAtt<-data.boot[,cov.ref]
covarRel<-data.boot[,cov.rel]
a <- as.matrix(covarAtt) %*% as.vector(as.numeric(coef.simul[b,]))

logVP.boot <- function(x)
{ 
return(sum(sapply(data.boot$indic[data.boot[,status.ref]==1], FUN = function(y) {
a[y] + as.matrix(covarRel[y,]) %*% as.vector(x[1:(dim(covarRel)[2])]) -
log(sum(exp(a[y:n]) * exp(as.matrix(covarRel[y:n,]) %*% as.vector(x[1:(dim(covarRel)[2])])))) } ) ) ) }    
    
model <- optim(init, logVP.boot, method = "Nelder-Mead", hessian = TRUE,
   control=list(fnscale=-1, maxit=100000))

        #iter<<-1
        #while(iter<=10) {temp.par<<-model$par; vrais<<-model$value
        #                 model<<-optim(temp.par, logVP.boot, method = "Nelder-Mead", hessian = TRUE,
        #                 control=list(fnscale=-1, maxit=100000))
        #                 iter<<-1*(model$value!=vrais) + (iter+1)*(model$value==vrais)
        #                 cat(model$value, "\n")}
        
        res.theta[b,]<-model$par
		#cat(b, "\n")
}

return(list(
matrix.coef = res.theta,
estim.coef = apply(res.theta, FUN="mean", MARGIN=2),
lower95.coef = apply(res.theta, FUN= function(x) {quantile(x, probs=0.025)}, MARGIN=2),
upper95.coef = apply(res.theta, FUN= function(x) {quantile(x, probs=0.975)}, MARGIN=2)
))

}