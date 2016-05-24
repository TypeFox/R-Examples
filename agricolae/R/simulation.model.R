`simulation.model` <-
function(model,file,categorical=NULL,k,console=FALSE) {
modelo<-model
parte<-strsplit(model,"~")[[1]]
model<-as.formula(model)
posicion<-0
if(length(categorical)>0){
posicion<-categorical
n<-length(posicion)
for( i in 1:n) {
pos<-posicion[i]
file[,pos]<-as.factor(file[,pos])
}
}
ecuacion<-lm(model,data=file)
xx<-data.frame(anova(ecuacion))
xx[,2]<-xx[,4]
fc<-xx[,4]
names(xx)<-c("Df","F value", "% Acceptance", "% Rejection", "Criterion")
gl<-anova(ecuacion)$Df
gk<-length(gl)-1
xx<-xx[-(gk+1),]
predicho<-predict(ecuacion)
m<-length(predicho)
sd.model<-sqrt(deviance(ecuacion)/gl[gk+1])
f<-rep(0,k)
cuenta <- rep(0,gk)
model <- paste("y","~",parte[2])
model<-as.formula(model)
for(i in 1:k){
errores<-rnorm(m,0,sd.model)
# Simula nuevos datos experimentales
y<-predicho+errores
simula<-lm(model,data=file)
for (j in 1:gk){
f[j]<-anova(simula)[j,4]
if(f[j] >= fc[j])cuenta[j]<-cuenta[j]+1
}
}
for( j in 1:gk){
xx[j,3]<-cuenta[j]*100/k
xx[j,4]<-100-xx[j,3]
if(xx[j,3]>50) xx[j,5]<- "acceptable"
if(xx[j,3]==50) xx[j,5]<- "uncertain"
if(xx[j,3]<50) xx[j,5]<- "nonacceptable"
}
if(console){
cat("\nSimulation of experiments\nUnder the normality assumption\n")
cat(rep("-",16),"\n")
cat("Proposed model:",modelo,"\n")
print(anova(ecuacion))
cat("---\n")
cat("Validation of the analysis of variancia for the proposed model\n")
cat("Simulations:",k,"\n\n")
print(xx)
cat("---\n\n")
}
out<-list(model=ecuacion, simulation=xx)
invisible(out)
}

