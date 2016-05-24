qLearn <-
function(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,s2Contrast,s1Contrast,alpha=0.05,bootNum=1000,s1Method="Regular",fixedXi,doubleBoot1Num=500,doubleBoot2Num=500,s1M,...){
s2_data<-subset(completeData,completeData[,s2Indicator]==1)
model<-getModel(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,...)
p<-model[[3]]
s2_var<-names(as.data.frame(model.matrix(s2Formula,data=s2_data,...)))[-1]
interaction<-paste(interact,s2Treat,sep=":")
interaction2<-paste(s2Treat,interact,sep=":")
reverse_match<-match(interaction2,s2_var)
if(sum(!is.na(reverse_match))>0){
match_num<-which(!is.na(reverse_match))
interaction[match_num]<-interaction2[match_num]
}
s2_var_noint<-s2_var[-match(s2Treat,s2_var)]
s2_var_noint<-s2_var_noint[-match(interaction,s2_var_noint)]
peudo_formula<- as.formula(paste(" ~ ", paste(s2_var_noint, collapse= "+")))
MM <- model.matrix(peudo_formula,data=completeData)
a<-match(names(as.data.frame(MM)),names(completeData))
if (sum(is.na(a))>0) {
MM <- MM[match(rownames(completeData),rownames(MM)),] 
rownames(MM) <- rownames(completeData)
extra<-MM[,is.na(a)]
completeData<-cbind(completeData,extra)
}
s1_var<-names(model.frame(s1Formula,data=completeData,...))

# stage 2 analysis
X2<-model.matrix(s2Formula,data=s2_data,...)
if (missing(s2Contrast)){
s2Contrast<-diag(ncol(X2))
}
stage2 <- matrix(0,nrow=dim(t(s2Contrast))[1],ncol=3)
s2_cf <- coef(model[[1]])
n2 <- dim(s2_data)[1]
bootest <- matrix(0,nrow=dim(t(s2Contrast))[1],ncol=bootNum)
for (i in 1:bootNum) {
index <- sample(1:n2,n2,replace=TRUE)
bootsamp <- s2_data[index,]
bootest[,i] <- s2Contrast%*%(2*s2_cf-coef(lm(s2Formula,data=bootsamp,...)))
}
for (i in 1:nrow(t(s2Contrast))) {
stage2[i,2:3] <- quantile(bootest[i,], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
}
stage2[,1]<-s2Contrast%*%s2_cf
colnames(stage2)<-c("S2_Estimator", "Lower", "Upper")

# stage 1 analysis
X1<-model.matrix(s1Formula,data=completeData,...)
if (missing(s1Contrast)){
s1Contrast<-diag(ncol(X1))
}
stage1 <- matrix(0,nrow=dim(t(s1Contrast))[1],ncol=3)
s1_cf <- coef(model[[2]])

k <- 0
n1 <- dim(completeData)[1]

# Use regular bootstrap for stage 1 as Default
if (s1Method=="Fixed Xi"){
s1M <- ceiling(n1^(1-p*(fixedXi/(1+fixedXi))))
}else if (s1Method=="Double Bootstrap"){
s1M <- chooseMDoubleBootstrap(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,alpha,doubleBoot1Num,doubleBoot2Num,...)
}else if (missing(s1M)){
s1M<-n1
}

cat("chosen value of m = ", s1M, "\n")

bootest <- matrix(0,nrow=dim(t(s1Contrast))[1],ncol=bootNum)
while (k<bootNum) {
index <- sample(1:n1,s1M,replace=TRUE)
bootsamp <- completeData[index,]
bootsamp_s2 <- subset(bootsamp,bootsamp[,s2Indicator]==1)
stage2cf <- coef(lm(s2Formula,data=bootsamp_s2,...))
bootsamp[,s1_var[1]]<-bootsamp[,s1_var[1]]+stage2cf[1]+as.matrix(bootsamp[,s2_var_noint])%*%stage2cf[s2_var_noint]+abs(stage2cf[s2Treat]+as.matrix(bootsamp[,interact])%*%stage2cf[interaction])

if (sum(is.na(bootsamp[,s1_var[1]]))<s1M){
k <- k+1
bootest[,k] <- s1Contrast%*%(2*s1_cf-coef(lm(s1Formula,data=bootsamp,...)))}
}

for (i in 1:nrow(t(s1Contrast))) {
stage1[i,2:3] <- quantile(bootest[i,], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
}
stage1[,1]<-s1Contrast%*%s1_cf
colnames(stage1)<-c("S1_Estimator", "Lower", "Upper")

print(round(stage2,2))
print(round(stage1,2))
object <- list(s1Coefficients = s1_cf, s1Inference = stage1, s2Coefficients = s2_cf, 
        s2Inference = stage2,s1Size=s1M) 
        object$call <- match.call()
        object
}

