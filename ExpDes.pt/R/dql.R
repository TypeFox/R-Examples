dql <-
function(trat, linha, coluna, resp, quali=TRUE, mcomp='tukey', sigT=0.05, sigF=0.05) {

Trat<-factor(trat)
Linha<-factor(linha)
Coluna<-factor(coluna)
anava<-aov(resp~Trat+Linha+Coluna)
tab<-summary(anava)

colnames(tab[[1]])<-c('GL','SQ','QM','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Tratamento','Linha','Coluna','Residuo','Total')
cv<-round(sqrt(tab[[1]][4,3])/mean(resp)*100, 2)
tab[[1]][5,3]=' '
cat('------------------------------------------------------------------------\nQuadro da analise de variancia
------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')


#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------\nTeste de normalidade dos residuos (Shapiro-Wilk)\n')
cat('p-valor: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('ATENCAO: a 5% de significancia, os residuos nao podem ser considerados normais!
------------------------------------------------------------------------\n')}
else{cat('De acordo com o teste de Shapiro-Wilk a 5% de significancia, os residuos podem ser considerados normais.
------------------------------------------------------------------------\n')}

if(tab[[1]][1,5]<sigF){

if(quali==TRUE) {

  if(mcomp=='tukey') tukey(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='duncan')duncan(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='lsd')   lsd(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='lsdb')  lsdb(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='sk')    scottknott(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='snk')   snk(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='ccboot')ccboot(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
  if(mcomp=='ccf')   ccF(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],sigT)
#  if(mcomp=="dnt")  {if(length(cont)==0) stop('Informe o nome do tratamento controle!')
#                     else dunnett(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],cont=cont,proc="dnt",alpha=sigT)}
#  if(mcomp=="sddnt"){if(length(cont)==0) stop('Informe o nome do tratamento controle!')
#                     else dunnett(resp,Trat,tab[[1]][4,1],tab[[1]][4,2],cont=cont,proc="sddnt",alpha=sigT)}
                     }
                 
else reg<-reg.poly(resp, trat, tab[[1]][4,1], tab[[1]][4,2], tab[[1]][1,1], tab[[1]][1,2])

                       }
else {
    cat('\nDe acordo com o teste F, as medias nao podem ser consideradas diferentes.\n')
mean.table<-tapply.stat(resp,trat,mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------\n')
}                       
residuals<-anava$residuals
coefficients<-anava$coefficients
effects<-anava$effects
fitted.values<-anava$fitted.values
means.trat<-tapply.stat(resp,trat,mean)
if(quali==FALSE) {
invisible(list(residuals=residuals, means.trat=means.trat, coefficients=coefficients, effects=effects, fitted.values=fitted.values,
'Regressao'=reg)) }

}
