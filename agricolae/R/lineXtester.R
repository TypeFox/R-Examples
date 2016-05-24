`lineXtester` <-
function(replications,lines,testers,y)
{
name.y <- deparse(substitute(y))

cat("\nANALYSIS LINE x TESTER: ", name.y, "\n")
replications<-as.factor(replications)
lines<-as.factor(lines)
testers<-as.factor(testers)
datos<-data.frame(Replications=replications,Lines=lines,Testers=testers,Y=y)
l<- length(levels(datos[,2]))
r<- length(levels(datos[,1]))
t<- length(levels(datos[,3]))
Treatments<-as.factor(paste(datos[,2],datos[,3]))
modelo1<-aov(Y ~ Replications + Treatments,data=datos)
matriz1<-as.matrix(anova(modelo1))
# line x tester
modelo4<- aov(Y~Lines*Testers,data=datos)
matriz4<-as.matrix(anova(modelo4))
# Crosses
datos2<-na.omit(datos)
mm<-tapply(datos2[,4], datos2[,2:3],mean,na.rm=TRUE)
cmm<-ncol(mm) ;rmm<-nrow(mm)
# Estimacion de SCA
SCA<-mm
for( i in 1:rmm) {
for( j in 1:cmm) {
SCA[i,j]<- round(mm[i,j]-mean(mm[,j],na.rm = TRUE)-mean(mm[i,],na.rm = TRUE)+mean(mm,na.rm = TRUE),3)
}
}
mm<-tapply(datos2[,4], datos2[,2],mean,na.rm=TRUE)
GCA.lines<-round(mm-mean(datos2[,4],na.rm=TRUE),3)
mm<-tapply(datos2[,4], datos2[,3],mean,na.rm=TRUE)
GCA.testers<-round(mm-mean(datos2[,4],na.rm=TRUE),3)
# ---------
Crosses<-as.factor(paste(datos2[,2],datos2[,3]))
modelo3<-aov(Y ~ Crosses,data=datos2)
matriz3<-as.matrix(anova(modelo3))
# Parents
datos3<-subset(datos,is.na(datos[,2]) | is.na(datos[,3]) )
Parents<-as.factor(paste(datos3[,2],datos3[,3]))
modelo2<-aov(Y ~ Parents,data=datos3)
matriz2<-as.matrix(anova(modelo2))
matriz5<-matriz1[2,]-matriz2[1,]-matriz3[1,]
matriz<-rbind(matriz1[1:2,],matriz2[1,],matriz5,matriz3[1,],matriz4[1:3,],matriz1[3,])
total1<-sum(matriz1[,1])
total2<-sum(matriz1[,2])
matriz5<-c(total1,total2,NA,NA,NA)
matriz<-rbind(matriz,matriz5)
for (i in 1:9){
matriz[i,3]<-matriz[i,2]/matriz[i,1]
matriz[i,4]<-round(matriz[i,3]/matriz[9,3],3)
matriz[i,5]<-round(1-pf(matriz[i,4],matriz[i,1],matriz[9,1]),4)
if (i==6 | i==7) {
matriz[i,4]<-round(matriz[i,3]/matriz[8,3],3)
matriz[i,5]<-round(1-pf(matriz[i,4],matriz[i,1],matriz[8,1]),4)
}
}
matriz[9,4]<-NA;matriz[9,5]<-NA
rownames(matriz)<-c("Replications","Treatments","Parents","Parents vs. Crosses","Crosses",
"Lines","Testers","Lines X Testers","Error","Total")
# Standard Errors for comnining ability effects
cm<-matriz[9,3]
s1<-sqrt(cm/(r*t))
s2<-sqrt(cm/(r*l))
s3<-sqrt(cm/r)
s4<-sqrt(2*cm/(r*t))
s5<-sqrt(2*cm/(r*l))
s6<-sqrt(2*cm/r)
# componentes geneticos
cov1<-(matriz[6,3]-matriz[8,3])/(r*t)
cov2<-(matriz[7,3]-matriz[8,3])/(r*l)
cov3<-(((l-1)*matriz[6,3]+(t-1)*matriz[7,3])/(l+t-2) - matriz[8,3])/(r*(2*l*t-l-t))
#      (((l-1)*2579.59    +(t-1)*859.463    )/(l+t-2) - 1770.296)   /(r*(2*l*t-l-t))
cov4<-((matriz[6,3]-matriz[9,3])+(matriz[7,3]-matriz[9,3])+(matriz[8,3]-matriz[9,3]))/(3*r)
+ (6*r*cov3-r*(l+t)*cov3)/(3*r)
# variancias
# Para F = 0 y F= 1
F<-0
var.A0 <-  cov3*(4/(1+F))
var.D0 <- ((matriz[8,3]-matriz[9,3])/r)*(2/(1+F))
F<-1
var.A1 <-  cov3*(4/(1+F))
var.D1 <- ((matriz[8,3]-matriz[9,3])/r)*(2/(1+F))
# Contribucion proporcional
c1<- matriz[6,2]*100/matriz[5,2]
c2<- matriz[7,2]*100/matriz[5,2]
c3<- matriz[8,2]*100/matriz[5,2]
# Impresion de resultados
cat("\nANOVA with parents and crosses",
    "\n==============================\n")
matriz1<-matriz[c(1,2,3,4,5,9,10),]
print(matriz1,na.print = "")
cat("\nANOVA for line X tester analysis",
    "\n================================\n")
matriz1<-matriz[6:9,]
print(matriz1,na.print = "")
cat("\nANOVA for line X tester analysis including parents",
    "\n==================================================\n")
print(matriz,na.print = "")
cat("\nGCA Effects:",
    "\n===========")
cat("\nLines Effects:\n")
print(GCA.lines)
cat("\nTesters Effects:\n")
print(GCA.testers)
cat("\nSCA Effects:",
    "\n===========\n")
print(SCA)
cat("\nStandard Errors for Combining Ability Effects:",
    "\n=============================================")
cat("\nS.E. (gca for line)   :",s1)
cat("\nS.E. (gca for tester) :",s2)
cat("\nS.E. (sca effect)     :",s3)
cat("\nS.E. (gi - gj)line    :",s4)
cat("\nS.E. (gi - gj)tester  :",s5)
cat("\nS.E. (sij - skl)tester:",s6,"\n")
cat("\nGenetic Components:",
    "\n==================")
cat("\nCov H.S. (line)   :",cov1)
cat("\nCov H.S. (tester) :",cov2)
cat("\nCov H.S. (average):",cov3)
cat("\nCov F.S. (average):",cov4)
cat("\nF = 0, Adittive genetic variance:",var.A0)
cat("\nF = 1, Adittive genetic variance:",var.A1)
cat("\nF = 0, Variance due to Dominance:",var.D0)
cat("\nF = 1, Variance due to Dominance:",var.D1,"\n")
cat("\nProportional contribution of lines, testers",
    "\n and their interactions to total variance",
    "\n===========================================")
cat("\nContributions of lines  :",c1)
cat("\nContributions of testers:",c2)
cat("\nContributions of lxt    :",c3,"\n")
return(modelo1)
}

