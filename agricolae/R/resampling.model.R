`resampling.model` <-
function(model,data,k,console=FALSE) {
modelo<-model
parte<-strsplit(model,"~")[[1]]
model<-as.formula(model)
ecuacion<-lm(model,data=data)
xx<-data.frame(anova(ecuacion),NA)
yy<-ecuacion$model[[1]]
fc<-xx[,4]
names(xx)<-c("Df","Sum Sq","Mean Sq","F value","Pr(>F)","Resampling")
m<-nrow(data)
gk<-nrow(xx)-1
f<-rep(0,gk)
cuenta <- rep(0,gk)
# Start Resampling
model <- paste("y","~",parte[2])
model<-as.formula(model)
for(i in 1:k){
y<-sample(yy,m)
resample<-lm(model,data=data)
for (j in 1:gk){
f[j]<-anova(resample)[j,4]
if(f[j] >= fc[j])cuenta[j]<-cuenta[j]+1
}
}
# finish resampling

for( j in 1:gk){
xx[j,6]<-cuenta[j]/k
}
if(console){
cat("\nResampling of the experiments\n")
cat(rep("-",14),"\n")
cat("Proposed model:",modelo,"\n")
cat("---\n")
cat("Resampling of the analysis of variancia for the proposed model\n")
cat("Determination of the P-Value by Resampling\n")
cat("Samples:",k,"\n\n")
xx<-as.matrix(xx)
print(xx,na.print="")
cat("---\n\n")
}
out<-list(model=resample, solution=xx,acum=cuenta,samples=k)
invisible(out)
}

