fat2.dbc <-
function(fator1, fator2, bloco, resp, quali=c(TRUE,TRUE), mcomp='tukey', fac.names=c('F1','F2'), sigT=0.05, sigF=0.05) {


cat('------------------------------------------------------------------------\nLegenda:\n')
cat('FATOR 1: ',fac.names[1],'\n')
cat('FATOR 2: ',fac.names[2],'\n------------------------------------------------------------------------\n\n')


fatores<-cbind(fator1,fator2)
Fator1<-factor(fator1)
Fator2<-factor(fator2)
Bloco<-factor(bloco)
nv1<-length(summary(Fator1))   #Diz quantos niveis tem o fator 1.
nv2<-length(summary(Fator2))   #Diz quantos niveis tem o fator 2.
J<-length(summary(Bloco))
lf1<-levels(Fator1)
lf2<-levels(Fator2)
lfB<-levels(Bloco)

anava<-aov(resp~ Bloco + Fator1*Fator2)
tab<-summary(anava)
#glB<-J-1
#glF1<-nv1-1
#glF2<-nv2-1
#glF1F2<-glF1*glF2
#glT<-J*nv1*nv2 -1
#glRes<-glT - glB - glF1 - glF2 - glF1F2 

#C<-sum(resp)^2/(J*nv1*nv2)
#SQB<-0
#for(i in 1:J) SQB<-SQB+sum(resp[Bloco==lfB[i]])^2
#SQB<-SQB/(nv1*nv2) - C






colnames(tab[[1]])<-c('GL','SQ','QM','Fc','Pr>Fc')
tab[[1]]<-rbind(tab[[1]],c(apply(tab[[1]],2,sum)))
rownames(tab[[1]])<-c('Bloco',fac.names[1],fac.names[2],paste(fac.names[1],'*',fac.names[2],sep=''),'Residuo','Total')
cv<-round(sqrt(tab[[1]][5,3])/mean(resp)*100, 2)
tab[[1]][6,3]=' '
cat('\nQuadro da analise de variancia\n------------------------------------------------------------------------\n')
print(tab[[1]])
cat('------------------------------------------------------------------------\nCV =',cv,'%\n')

#Teste de normalidade
pvalor.shapiro<-shapiro.test(anava$residuals)$p.value
cat('\n------------------------------------------------------------------------
Teste de normalidade dos residuos (Shapiro-Wilk)\n')
cat('p-valor: ',pvalor.shapiro, '\n')
if(pvalor.shapiro<0.05){cat('ATENCAO: a 5% de significancia, os residuos nao podem ser considerados normais!
------------------------------------------------------------------------\n')}
else{cat('De acordo com o teste de Shapiro-Wilk a 5% de significancia, os residuos podem ser considerados normais.
------------------------------------------------------------------------\n')}

#Para interacao nao significativa, fazer...
if(tab[[1]][4,5]>sigF) {                                    
cat('\nInteracao nao significativa: analisando os efeitos simples
------------------------------------------------------------------------\n')
fatores<-data.frame('fator 1'=fator1,'fator 2' = fator2)

for(i in 1:2){

#Para os fatores QUALITATIVOS, teste de Tukey
if(quali[i]==TRUE && tab[[1]][i+1,5]<=sigF) {                  
    cat(fac.names[i])
      if(mcomp=='tukey'){
    tukey(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='duncan'){
    duncan(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='lsd'){
    lsd(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='lsdb'){
    lsdb(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='sk'){
    scottknott(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=='snk'){
    snk(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                    }
  if(mcomp=="ccboot"){
  ccboot(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                     }
  if(mcomp=="ccf"){
  ccF(resp,fatores[,i],tab[[1]][5,1],tab[[1]][5,2],sigT)
                     }
                   }
if(quali[i]==TRUE && tab[[1]][i+1,5]>sigF) {
    cat(fac.names[i])
    cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
}

#Para os fatores QUANTITATIVOS, regressao
if(quali[i]==FALSE && tab[[1]][i+1,5]<=sigF){
    cat(fac.names[i])
    reg.poly(resp, fatores[,i], tab[[1]][5,1], tab[[1]][5,2], tab[[1]][i+1,1], tab[[1]][i+1,2])
}

if(quali[i]==FALSE && tab[[1]][i+1,5]>sigF) {
    cat(fac.names[i])
    cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n\n')
mean.table<-tapply.stat(resp,fatores[,i],mean)
colnames(mean.table)<-c('Niveis','Medias')
print(mean.table)
cat('------------------------------------------------------------------------')
                            }

cat('\n')
}

}

#Se a interacao for significativa, desdobrar a interacao
if(tab[[1]][4,5]<=sigF){
cat("\n\n\nInteracao significativa: desdobrando a interacao
------------------------------------------------------------------------\n")

#Desdobramento de FATOR 1 dentro do niveis de FATOR 2
cat("\nDesdobrando ", fac.names[1], ' dentro de cada nivel de ', fac.names[2], '
------------------------------------------------------------------------\n')

des1<-aov(resp~ Bloco + Fator2/Fator1)

l1<-vector('list',nv2)
names(l1)<-names(summary(Fator2))
v<-numeric(0)
for(j in 1:nv2) {
        for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
        l1[[j]]<-v
        v<-numeric(0)
                }
des1.tab<-summary(des1,split=list('Fator2:Fator1'=l1))[[1]]
#Montando a tabela de ANAVA do des1
glB=tab[[1]][1,1]
glb=nv2-1
glf1=c(as.numeric(des1.tab[4:(nv2+3),1]))
glE=tab[[1]][5,1]
glT=tab[[1]][6,1]

SQB=tab[[1]][1,2]
SQb=tab[[1]][3,2]
SQf1=c(as.numeric(des1.tab[4:(nv2+3),2]))
SQE=tab[[1]][5,2]
SQT=tab[[1]][6,2]

QMB=SQB/glB
QMb=SQb/glb
QMf1=SQf1/glf1
QME=SQE/glE
QMT=SQT/glT

FcB=QMB/QME
Fcb=QMb/QME
Fcf1=QMf1/QME

rn<-numeric(0)
for(j in 1:nv2){ rn<-c(rn, paste(paste(fac.names[2],':',fac.names[1],sep=''),lf2[j]))}

anavad1<-data.frame("GL"=c(round(c(glB, glb, glf1, glE, glT))),
"SQ"=c(round(c(SQB,SQb,SQf1,SQE,SQT),5)),
"QM"=c(round(c(QMB,QMb,QMf1,QME,QMT),5)),
"Fc"=c(round(c(FcB,Fcb,Fcf1),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fcb,glb,glE),1-pf(Fcf1,glf1,glE)),4),' ', ' '))
rownames(anavad1)=c("Bloco",fac.names[2],rn,"Residuo","Total")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad1)
cat('------------------------------------------------------------------------\n\n')


for(i in 1:nv2) {
  if(des1.tab[(i+3),5]<=sigF){
    if(quali[1]==TRUE){
                      cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
                        if(mcomp=='tukey'){
                          tukey(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                          }
                        if(mcomp=='duncan'){
                          duncan(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=='lsd'){
                          lsd(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=="ccboot"){
                        ccboot(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=="ccF"){
                        ccF(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                      }
    else{  #regressao
    cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'
------------------------------------------------------------------------')
    reg.poly(resp[Fator2==lf2[i]], fator1[Fator2==lf2[i]], tab[[1]][5,1], tab[[1]][5,2], des1.tab[i+3,1], des1.tab[i+3,2][[1]])
        }
                              }
    else{cat('\n\n',fac.names[1],' dentro do nivel ',lf2[i],' de ',fac.names[2],'\n')
    cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
        mean.table<-tapply.stat(resp[Fator2==lf2[i]],fatores[,1][Fator2==lf2[i]],mean)
        colnames(mean.table)<-c('  Niveis','    Medias')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }
                 }
cat('\n\n')

#Desdobramento de FATOR 2 dentro do niveis de FATOR 1
cat("\nDesdobrando ", fac.names[2], ' dentro de cada nivel de ', fac.names[1], '
------------------------------------------------------------------------\n')

des2<-aov(resp~ Bloco + Fator1/Fator2)

l2<-vector('list',nv1)
names(l2)<-names(summary(Fator1))
v<-numeric(0)
for(j in 1:nv1) {
        for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
        l2[[j]]<-v
        v<-numeric(0)
                }
des2.tab<-summary(des2,split=list('Fator1:Fator2'=l2))[[1]]
#Montando a tabela de ANAVA do des2
gla=nv1-1
glf2=c(as.numeric(des2.tab[4:(nv1+3),1]))

SQa=tab[[1]][2,2]
SQf2=c(as.numeric(des2.tab[4:(nv1+3),2]))

QMa=SQa/gla
QMf2=SQf2/glf2

Fca=QMa/QME
Fcf2=QMf2/QME

rn<-numeric(0)
for(i in 1:nv1){ rn<-c(rn, paste(paste(fac.names[1],':',fac.names[2],sep=''),lf1[i]))}

anavad2<-data.frame("GL"=c(round(c(glB, gla, glf2, glE, glT))),
"SQ"=c(round(c(SQB,SQa,SQf2,SQE,SQT),5)),
"QM"=c(round(c(QMB,QMa,QMf2,QME,QMT),5)),
"Fc"=c(round(c(FcB,Fca,Fcf2),4),'',''),
"Pr>Fc"=c(round(c(1-pf(FcB,glB,glE),1-pf(Fca,gla,glE),1-pf(Fcf2,glf2,glE)),4),' ', ' '))
rownames(anavad2)=c("Bloco",fac.names[1],rn,"Residuo","Total")
cat('------------------------------------------------------------------------
Quadro da analise de variancia\n------------------------------------------------------------------------\n')
print(anavad2)
cat('------------------------------------------------------------------------\n\n')


for(i in 1:nv1) {
  if(des2.tab[(i+3),5]<=sigF){
    if(quali[2]==TRUE){
                      cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
                      if(mcomp=='tukey'){
                          tukey(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                          }
                        if(mcomp=='duncan'){
                          duncan(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=='lsd'){
                          lsd(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=='lsdb'){
                          lsdb(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                         }
                        if(mcomp=='sk'){
                          scottknott(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                       }
                        if(mcomp=='snk'){
                          snk(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                        }
                        if(mcomp=="ccboot"){
                        ccboot(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                        if(mcomp=="ccf"){
                        ccF(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],tab[[1]][5,1],tab[[1]][5,2],sigT)
                                           }
                      }
    else{  #regressao
        cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'
------------------------------------------------------------------------')
        reg.poly(resp[Fator1==lf1[i]], fator2[Fator1==lf1[i]], tab[[1]][5,1], tab[[1]][5,2], des2.tab[i+3,1], des2.tab[i+3,2][[1]])
        }
                             }
    else{cat('\n\n',fac.names[2],' dentro do nivel ',lf1[i],' de ',fac.names[1],'\n')
    cat('\nDe acordo com o teste F, as medias desse fator sao estatisticamente iguais.\n')
        mean.table<-tapply.stat(resp[Fator1==lf1[i]],fatores[,2][Fator1==lf1[i]],mean)
        colnames(mean.table)<-c('  Niveis','    Medias')
        print(mean.table)
        cat('------------------------------------------------------------------------\n')
        }
                              
                }
}


}
