reg.poly <-
function(resp, trat, glres, SQres, gltrat, SQtrat) { 

QMres<-SQres/glres

cat('\nAjuste de modelos polinomiais de regressao\n------------------------------------------------------------------------\n')

X<-matrix(1,length(trat),4)
X[,2]<-trat
X[,3]<-trat^2
X[,4]<-trat^3

mean.table<-tapply.stat(resp,trat,mean)
colnames(mean.table)<-c('  Niveis','   Medias Observadas')

###############################################################################
#reta
###############################################################################
b=ginv(t(X[,1:2])%*%X[,1:2], tol=.Machine$double.eps)%*%t(X[,1:2])%*%resp
ep = sqrt(diag(ginv(t(X[,1:2])%*%X[,1:2], tol=.Machine$double.eps)*QMres))
tc = b/ep
pv = 2*pt(abs(tc),glres,lower.tail=FALSE)     
tm1<-data.frame('Estimativa' = round(b,8),'Erro padrao' = round(ep,5),'tc'=round(tc,5),'p-valor' = round(pv,5))
rownames(tm1)<-c('b0','b1')

aov.m1<-anova(lm(resp~trat))
if(dim(mean.table)[1]==2){r2m1<-1}
if(dim(mean.table)[1]>2) {r2m1<-aov.m1[1,2]/SQtrat}

#ANAVA da regressao linear
nomes1<-c("Efeito linear","Desvios de Regressao","Residuos")
anava1<-data.frame('GL'=c(1,(gld=c(gltrat-1)),(glr=glres)),
                  'SQ'=c(round(c(aov.m1[[2]][1],(sqd=c(SQtrat-aov.m1[[2]][1])),SQres),5)),
                  'QM'=c(round(c(aov.m1[[3]][1],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), QMres),5)),
                  'Fc'=c(round(c((fcl=aov.m1[[3]][1]/QMres),(fc=qmd/QMres)),2),''),
                  'p-valor'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1}else{pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava1)<-nomes1

output1<-list('Modelo linear
------------------------------------------------------------------------' = tm1,
              'R2 do modelo linear' = r2m1,
              'Analise de variancia do modelo linear' = anava1)
print(output1,right=TRUE)
cat('------------------------------------------------------------------------\n')
###############################################################################
  if(dim(mean.table)[1]>2) {
###############################################################################
#parabola
###############################################################################
b2=ginv(t(X[,1:3])%*%X[,1:3], tol=.Machine$double.eps)%*%t(X[,1:3])%*%resp
ep2 = sqrt(diag(ginv(t(X[,1:3])%*%X[,1:3], tol=.Machine$double.eps)*QMres))
tc2 = b2/ep2
pv2 = 2*pt(abs(tc2),glres,lower.tail=FALSE)  
tm2<-data.frame('Estimativa' = round(b2,8),'Erro padrao' = round(ep2,5),'tc'=round(tc2,5),'p-valor' = round(pv2,5))
rownames(tm2)<-c('b0','b1','b2')

t2<-trat^2
aov.m2<-anova(lm(resp~trat+t2))
if(dim(mean.table)[1]==3){r2m2<-1}
if(dim(mean.table)[1]>3) {r2m2<-(aov.m2[1,2]+aov.m2[2,2])/SQtrat}

#ANAVA da regressao quadratica
nomes2<-c("Efeito linear","Efeito quadratico","Desvios de Regressao","Residuos")
anava2<-data.frame('GL'=c(aov.m2[[1]][1:2],(gld=c(gltrat-2)),(glr=glres)),
                  'SQ'=c(round(c(aov.m2[[2]][1:2],(sqd=c(SQtrat-sum(aov.m2[[2]][1:2]))),SQres),5)),
                  'QM'=c(round(c(aov.m2[[3]][1:2],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), QMres),5)),
                  'Fc'=c(round(c((fcl=aov.m2[[3]][1:2]/QMres),(fc=qmd/QMres)),2),''),
                  'p-valor'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1}else{pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava2)<-nomes2

output2<-list('Modelo quadratico
------------------------------------------------------------------------' = tm2,
              'R2 do modelo quadratico' = r2m2,
              'Analise de variancia do modelo quadratico' = anava2)
print(output2,right=TRUE)
cat('------------------------------------------------------------------------\n')

                            }
###############################################################################
  if(dim(mean.table)[1]>3) {
###############################################################################
#cubica
###############################################################################
b3=ginv(t(X[,1:4])%*%X[,1:4], tol=.Machine$double.eps)%*%t(X[,1:4])%*%resp
ep3 = sqrt(diag(ginv(t(X[,1:4])%*%X[,1:4], tol=.Machine$double.eps)*QMres))
tc3 = b3/ep3
pv3 = 2*pt(abs(tc3),glres,lower.tail=FALSE)  
tm3<-data.frame('Estimativa' = round(b3,8),'Erro padrao' = round(ep3,5),'tc'=round(tc3,5),'p-valor' = round(pv3,5))
rownames(tm3)<-c('b0','b1','b2','b3')

t3<-trat^3
aov.m3<-anova(lm(resp~trat+t2+t3))
if(dim(mean.table)[1]==4){r2m3<-1}
if(dim(mean.table)[1]>4) {r2m3<-(aov.m3[1,2]+aov.m3[2,2]+aov.m3[3,2])/SQtrat}

#ANAVA da regressao cubica
nomes3<-c("Efeito linear","Efeito quadratico","Efeito cubico","Desvios de Regressao","Residuos")
anava3<-data.frame('GL'=c(aov.m3[[1]][1:3],(gld=c(gltrat-3)),(glr=glres)),
                  'SQ'=c(round(c(aov.m3[[2]][1:3],(sqd=c(SQtrat-sum(aov.m3[[2]][1:3]))),SQres),5)),
                  'QM'=c(round(c(aov.m3[[3]][1:3],(if(gld==0){qmd=0} else{qmd=(sqd/gld)}), QMres),5)),
                  'Fc'=c(round(c((fcl=aov.m3[[3]][1:3]/QMres),(fc=qmd/QMres)),2),''),
                  'p-valor'=c(round(c(pf(fcl,1,glr,lower.tail=FALSE),(if(gld==0){pv=1} else {pv=1-pf(fc,gld,glr)})),5),''))
rownames(anava3)<-nomes3

output3<-list('Modelo cubico
------------------------------------------------------------------------' = tm3,
             'R2 do modelo cubico' = r2m3, 
             'Analise de variancia do modelo cubico' = anava3
              )
print(output3,right=TRUE)
cat('------------------------------------------------------------------------\n')
                            }
###############################################################################
if(dim(mean.table)[1]>3) {return(list("Quadro de medias" = mean.table, "Coeficientes reta" = b, "R2" = r2m1, "Coeficientes parabola" = b2,
                                       "R2" = r2m2, "Coeficientes cubica" = b3, "R2" = r2m3)) }
if(dim(mean.table)[1]==3){return(list("Quadro de medias" = mean.table, "Coeficientes reta" = b, "R2" = r2m1, "Coeficientes parabola" = b2,
                                       "R2" = r2m2)) }
if(dim(mean.table)[1]<3) {return(list("Quadro de medias" = mean.table, "Coeficientes reta" = b, "R2" = r2m1)) }

cat('------------------------------------------------------------------------\n\n')

}
